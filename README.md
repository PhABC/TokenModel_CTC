# TokenModel_CTC
Token task model for neural data by Cisek, Thivierge and Castonguay
----------------------------------------------------------------------------------------


This model explores how the decision manifold is shaped by different neural architecture.
It is also based on existing data (from PMd and M1 in monkeys) that should be published
by the end of 2016.


# HOW TO USE ----------------------------------------------------------------------------

In order to play with the different architecture, there are some parameters more imp
than others. Here, it will be describe which parameters are the most relevant and
what they mean.

# Seed
Seed allows you to replay a simulation, which can be useful if testing parameters.
	+ seed  = Saving the current seed. Comment out if you want to use the previous seed

# Stimulation parameters
In this section, you can choose general parameters for simulation:

	+ N     = Number of neurons.
	+ T     = Time of simulation.
	+ onset = Stimulus onset.

# Connection parametes
These parameters are the key parameters to determine the architecture of the network.
There are connection parameters from between regions and within regions.The inhibition
and excitation strenght are proportional since the weight matrix is only a single sunken
gaussian curve. The standart deviation, the sunking ratio and the amplitude can be
choosen with the following parameters:

	+ Ww_*%   = Magnitude of the matrix from region * to region %.
	+ Sunk_*% = Proporting of excitory to inhibitory units from * to %.
		If set to 1, then all the weights are negative.
		If set to 0, all the weights are positive.
	+ Wsd_*%  = Standart deviation of the weight matrix from * to %.

#Input parameters
In this section, you can select what will be fed to the neurons:

	Stimuli parameters
	+ c     = type of stimuli with the parameter.
		1 : easy trial, 2 : misleading trial, 3 : ambiguous trials
	+ stimW = The amplitude of the stimuli.

	Biais
	+ bias = Is a biais input to share by all neurons.

	Noise
	+ fg = Fast gaussian noise strength, unique to each neurons.
	+ sg = Slow gaussian noise strength, share with all neurons.

	Urgency
	+ Utype = The type of urgency signal. 1 = additive urgency & 2 = multiplicative urgency
	+ Uw    = The amplitude of the urgency signal
	+ Uslop = Slope of the urgency signal (in average)


