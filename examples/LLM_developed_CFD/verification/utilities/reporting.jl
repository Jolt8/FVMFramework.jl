struct VerificationResult
    name::String
    component::String
    status::Symbol
    summary::String
    metrics::Dict{String, Any}
    expected::Dict{String, Any}
    diagnostics::Dict{String, Any}
    elapsed_seconds::Float64
end

struct VerificationReport
    target::String
    level::Symbol
    started_at::DateTime
    finished_at::DateTime
    environment::Dict{String, Any}
    results::Vector{VerificationResult}
end

function _verification_result(
    name,
    component,
    passed,
    summary;
    metrics = Dict{String, Any}(),
    expected = Dict{String, Any}(),
    diagnostics = Dict{String, Any}(),
    elapsed_seconds = 0.0,
)
    if passed
        status = :pass
    else
        status = :fail
    end
    return VerificationResult(
        name,
        component,
        status,
        summary,
        metrics,
        expected,
        diagnostics,
        elapsed_seconds,
    )
end

function _run_check(check, name, component)
    start_time = time()
    try
        result = check()
        return _verification_result(
            name,
            component,
            result.passed,
            result.summary;
            metrics = get(result, :metrics, Dict{String, Any}()),
            expected = get(result, :expected, Dict{String, Any}()),
            diagnostics = get(result, :diagnostics, Dict{String, Any}()),
            elapsed_seconds = time() - start_time,
        )
    catch exception
        diagnostic = sprint(showerror, exception, catch_backtrace())
        return VerificationResult(
            name,
            component,
            :error,
            "check raised an exception",
            Dict{String, Any}(),
            Dict{String, Any}(),
            Dict{String, Any}(
                "exception_type" => string(typeof(exception)),
                "exception" => diagnostic,
            ),
            time() - start_time,
        )
    end
end

function has_failures(report::VerificationReport)
    return any(result -> result.status != :pass, report.results)
end

function _result_dict(result::VerificationResult)
    return Dict{String, Any}(
        "name" => result.name,
        "component" => result.component,
        "status" => string(result.status),
        "summary" => result.summary,
        "metrics" => result.metrics,
        "expected" => result.expected,
        "diagnostics" => result.diagnostics,
        "elapsed_seconds" => result.elapsed_seconds,
    )
end

function report_dict(report::VerificationReport)
    pass_count = count(result -> result.status == :pass, report.results)
    fail_count = count(result -> result.status == :fail, report.results)
    error_count = count(result -> result.status == :error, report.results)
    if has_failures(report)
        overall_status = "fail"
    else
        overall_status = "pass"
    end

    return Dict{String, Any}(
        "schema_version" => "1.0",
        "target" => report.target,
        "level" => string(report.level),
        "started_at" => string(report.started_at),
        "finished_at" => string(report.finished_at),
        "status" => overall_status,
        "summary" => Dict{String, Any}(
            "passed" => pass_count,
            "failed" => fail_count,
            "errors" => error_count,
            "total" => length(report.results),
        ),
        "environment" => report.environment,
        "results" => [_result_dict(result) for result in report.results],
    )
end

function _json_escape(value)
    return replace(
        string(value),
        '\\' => "\\\\",
        '"' => "\\\"",
        '\n' => "\\n",
        '\r' => "\\r",
        '\t' => "\\t",
    )
end

function _write_json(io, value, indent, depth)
    indentation = repeat(" ", indent * depth)
    child_indentation = repeat(" ", indent * (depth + 1))

    if value === nothing || value === missing
        print(io, "null")
    elseif value isa Bool
        if value
            print(io, "true")
        else
            print(io, "false")
        end
    elseif value isa Integer
        print(io, value)
    elseif value isa AbstractFloat
        if isfinite(value)
            print(io, value)
        else
            print(io, '"', value, '"')
        end
    elseif value isa AbstractString || value isa Symbol || value isa DateTime || value isa VersionNumber
        print(io, '"', _json_escape(value), '"')
    elseif value isa NamedTuple
        _write_json(io, Dict(string(key) => item for (key, item) in pairs(value)), indent, depth)
    elseif value isa AbstractDict
        entries = sort!(collect(pairs(value)); by = pair -> string(first(pair)))
        print(io, '{')
        if !isempty(entries)
            print(io, '\n')
            for (entry_index, (key, item)) in enumerate(entries)
                print(io, child_indentation, '"', _json_escape(key), "\": ")
                _write_json(io, item, indent, depth + 1)
                if entry_index < length(entries)
                    print(io, ',')
                end
                print(io, '\n')
            end
            print(io, indentation)
        end
        print(io, '}')
    elseif value isa Tuple || value isa AbstractArray
        print(io, '[')
        if !isempty(value)
            print(io, '\n')
            for (item_index, item) in enumerate(value)
                print(io, child_indentation)
                _write_json(io, item, indent, depth + 1)
                if item_index < length(value)
                    print(io, ',')
                end
                print(io, '\n')
            end
            print(io, indentation)
        end
        print(io, ']')
    else
        print(io, '"', _json_escape(value), '"')
    end
end

function report_json(report::VerificationReport; indent = 2)
    io = IOBuffer()
    _write_json(io, report_dict(report), indent, 0)
    return String(take!(io))
end

function _short_diagnostic(value)
    rendered = sprint(show, value)
    if length(rendered) > 300
        return first(rendered, 297) * "..."
    end
    return rendered
end

function report_string(report::VerificationReport)
    io = IOBuffer()
    println(io, "FVMFramework Verification Report")
    println(io, "Target: $(report.target)    Level: $(report.level)")

    components = unique(result.component for result in report.results)
    for component in components
        println(io)
        println(io, uppercase(component))
        for result in report.results
            if result.component != component
                continue
            end
            status = uppercase(string(result.status))
            @printf(io, "%-5s %-34s %s\n", status, result.name, result.summary)
            if result.status != :pass
                if !isempty(result.metrics)
                    println(io, "      metrics: $(_short_diagnostic(result.metrics))")
                end
                if !isempty(result.expected)
                    println(io, "      expected: $(_short_diagnostic(result.expected))")
                end
                for (key, value) in sort!(collect(result.diagnostics); by = first)
                    println(io, "      $key: $(_short_diagnostic(value))")
                end
            end
        end
    end

    pass_count = count(result -> result.status == :pass, report.results)
    if has_failures(report)
        overall_status = "FAIL"
    else
        overall_status = "PASS"
    end
    println(io)
    println(io, "RESULT: $overall_status ($pass_count/$(length(report.results)) checks passed)")
    return String(take!(io))
end

function write_report(report::VerificationReport, output_directory)
    mkpath(output_directory)
    text_path = joinpath(output_directory, "verification_report.txt")
    json_path = joinpath(output_directory, "verification_report.json")
    write(text_path, report_string(report))
    write(json_path, report_json(report))
    return (text = text_path, json = json_path)
end
