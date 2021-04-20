# Go Corona Go : Nonlinear Dynamics based Study to build an Interventional Model
Proposing a Network Centrality based Intervention Model on basis of understanding various epidemic models, starting from classical SIR model with / without vital dynamics to models under the influence of intervention, saturated infection, quadratic infection and vaccination. 
### Study by : Sayantan Jana and Sarthak Singhal 

#### Abstract
Throughout history, infectious diseases, which are illnesses caused by a disease agent that can be
transmitted from organism to organism, have had a big impact on the human health and social
development through their epidemic course of action. Presently the world is under the COVID
pandemic which started as an outbreak in China. This causes large scale harm to the society and
therefore it is imperative to study such diseases to prevent them from growing to a stage where it
is uncontrollable. Thus these diseases are modelled and we study their evolution in a hope to create
some intervention by either vaccinating or lockdown or other means. Through this work we aim to
study an epidemic model with the effects of various forces which may intensify/suppress
the epidemic, mostly following the work of Wendi Wang and try to build an intervention
model on top of the prevention models at local spaces, using percolation centrality as the
basis .

### Brief Work Description 
- We started with understanding the very basics of SIR models, without and with vital dynamics,
understand what do SIR, SIRS, SEIR, SEIQR and SEIRS signify.
- Next we studied the survey paper ”Epidemic Models with Nonlinear Infection Forces” by Wendi
Weng, which had followed the previous work of Ruan and Weng, where we learnt of the effect
of forces concerned with the intervention as well as a saturated infection force, studied the
proofs of existence of the equilibrium points, disease free equilibrium and endemic equilibrium.
- Constructing time course simulations using Xppaut for the coupled system of nonlinear ODEs
in each of the above models, we observed their dynamics.
- We studied and discussed briefly about percolation centrality, as proposed by Piraveenan,
Prokopenko, Hossain and using that and the nonlinear dynamic models we studied so far,
we propose an intervention model based on percolation centrality algorithm for a large area
consisting of multiple smaller zones. However, the computation would require O(N^3 · S), where
S is constant time to run the simulations on NLD model. Hence, we have implemented a parallel
percolation centrality computation algorithm in CUDA.
- During the course of reading for the project, we explored some topics which we didn’t go into
greater depth and hence have put in the future scope section.

The visualisations of the functions under study are put in the folder function_visualisations. All the implementation and Analysis visualisations have been put in code folder. For compiling,running and formatting for the individual codes, README has been provided separately for the folders inside code folder.

Note : We initially planned to just have focus on study of intervention models and later found that we could connected the graph theoretic centrality based techniques to this. 



