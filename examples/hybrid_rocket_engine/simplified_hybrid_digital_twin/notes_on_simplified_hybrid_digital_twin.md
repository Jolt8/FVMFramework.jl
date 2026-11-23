Overall, compared to the previous models that I built, I was certainly the most happy with this one. 

# Things I liked about this project

# An LLM realized that the injector could be merged to no longer require the mid section
- This is what I actually absolutely love LLMs for, they are extremely creative when it comes to proposing solutions that really wouldn't of ever occured to me
- Because we did this, we no longer had to deal with the extremely sensitive mid_section transient, making it significantly easier to solve the problem 

# No moving mesh was used
- Since this script did not use a moving mesh, it was significantly easier to get a working model and I never really had to rely on an LLM to code anything for me, allowing me to take full ownership of the model

# NASA CEA was used
- Compared to Arrhenius chemistry which was extremely stiff for this problem, NASA CEA was extremely practical here and worked very well

# Problems and Imporvements

# Using rocketcea and being forced to use HP mode
Problem:
- Since we're using rocketcea and it's forced to be in HP mode, it does not map especially cleanly to FVM

Improvements:
- Instead of using rocketcea in the future, we should ideally create our own CEA engine in Julia that utilizes its own gibbs or helmholtz free energy minimization routine combined with NASA's 9 coefficient property polynomials 
- Another option would just be to create an application that takes in the data from NASA CEARUN, but I would like to develop my own to make it significantly more customizable and so I understand it better