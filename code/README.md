Code for a project in the "Optimal Control" exam.

The theoretical model folder contains code for the simulation and control of an SIR model using the forwards backwards method as optimization solver. It also features a script to perform parameter identification from the fake measurement of a solution of the model to test its identification.

Since the identification performs very poorly with the chosen model another folder containing the same procedure for a simpler model is available.

Using this simpler model, real world data coming from the 2016 Ebola epidemic was analyzed to determine the parameters on the model. Finally, a model with these estimated parameters was controlled using the same method as above. 
