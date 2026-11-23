# Why was axisymmetric FVM the wrong approach this for this project?

Here are a couple issues that I remember from working on this simulation script:

# Problems and Improvements in the future:

# A moving fuel grain was used
Problem:
- Due to my ideal that I would one day be able to link the regression observed in the fuel grain to what happened during the test fire, I constantly wanted to use a model that could explicitly track this moving fuel grain over time. 
- In hindsight, this really was not the right choice for a few reasons:
    1. It adds an insane amount of extra complexity that made me rely on an LLM to get it working properly
    2. While I haven't confirmed this for sure, I don't think that the moving fuel grain would've had a significant impact on the information that we were trying to extract out of the model
- Since I used an LLM to get this implemented faster, I really lost ownership of certain parts of this model and I lost interest in it

Improvements:
- Don't preemptively build a model that can interface iwth future collected data because it's almost always better to do something when it's actually needed rather than because you think it might be useful in the future
- While it's hard to do this, always critically evaluate whehther or not a mechanism that is complicated actually adds that much to the problem and don't get blinded by the belief that (more accurate physics == always better)
- If you ever feel tempted to use an LLM for a systems engineering model, just DON'T, it's not worth it, I don't know if it's just because I'm overly petty and I hate when I lose ownership of my code, but it genuinely feels that for systems models even a slight amount of lost information makes everything significantly harder to work on 
- Perhaps we should develop a more generally applicable moving mesh methodology (even if it's generated with the help of an LLM) so that, even though it's abstracted, I get to choose where, when, and why something is or isn't used 
    - This might be impossible because moving meshes typically require that momentum and other conservation of energy implications be dealt with (at least for SLAU2)

# Arrenhius chemistry was used (not going to write this out fully, just see my explanation for the FVM script used in this project)