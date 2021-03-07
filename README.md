# Algorithms for Adaptive Influence Maximization

Codes for the EptAIM, i.e., the adaptive algorithm for adaptive IM problem that provides expected approximation ratio. More details refer to paper: https://arxiv.org/abs/2004.06469.



### Running Environment

Linux-based OS



### Input requirements

**Dataset files**: dataset files are referred to $\textit{Tested-Dataset}$ repository.

**Realizations files**: please refer to $\textit{Tools}$ repository for source code to generate realizations.

*Remark*: please note that generated realization files are supposed to locate in the same folder of dataset file according to my implementation.



### How to run

#### Running command
./algo -dataset path_to_dataset -model IC -epsilon ε -k seed_number -batch batch_size -seedfile filename -time time_number

#### Explain
--epsilon:  an float number in range (0,1) to control the approximation error.

--k: the number of seed nodes to be selected.

--batch: the size of batch b selected each time. 

--seedfile: the file records the k seed nodes selected.

--time: the number of the algorithm repeated.

#### For example:

./exp_epic -dataset dataset/hep/ -model IC -epsilon 0.5 -k 500 -batch 50 -seedfile seed -time 1



### Datasets

Tested datasets can be downloaded in $\textit{Tested-Dataset}$ repository.



### Remark

If there are any problems, please contact khuang005@ntu.edu.sg / kkhuang@nus.edu.sg.

