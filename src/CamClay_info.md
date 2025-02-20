## Information on the Cam Clay Model


These are some notes on the cam clay model at the importance of the equations that are in the model and how they should be formulated together


* Function for the yield surface

* Function for the normal to the yield surface

* Associative flow so the plastic potential is the yield surface

* Information on how the preconsolidation pressure is updated. Need to make sure the sign conventions work out. I think there is a problem with the sign convention or how the operations goes about assigning the updates and the minimums that are required

* Need to think about the objective stress rate and the implications it has for the MPM model. I don't want to use a specific stress rate if I don't have to. It would be preferable to rotate the materials back to their initial configuration and run the constitutive model on that

p_c: preconsolidation pressure
p: mean_stress, defined that mean stress is always positive
F: Value of the yield function
lambda: Slope of the normal compression line in semi-log space
kappa: Slope of the swelling line in semi-log space

* Need to be aware of the concerns about the mean stress and how that is going to impact the form of the yield surface and the subsequent derivatives

* Non-linear elasticity makes the model thermodynamically inconsistent because the elastic portion is stress path dependent

## Questions

* What impact does using different invariants have on calculating the derivative othe yield surface and the plastic potential? 
    - They are all related so it should just be an application of the chain rule but what implication does that have on the results you get?

* How can different methods for calculating the yield surface be handled?
    - Using classes would probably be the easiest way because that could be set.
        + This would have an impact on the computational speed of the model
        
