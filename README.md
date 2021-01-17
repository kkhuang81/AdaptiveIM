# Algorithms for Adaptive Influence Maximization

Codes for the EptAIM, i.e., the adaptive algorithm for adaptive IM problem that provides expected approximation ratio. More details refer to [paper]:https://arxiv.org/abs/2004.06469.

# Running Environment

Linux-based OS

# How to run

./algo -dataset path_to_dataset -model IC -epsilon $\epsilon$ -k seed_number -batch batch_size -seedfile filename -time time_number

## Explain
--epsilon:  an float number in range (0,1) to control the approximation error.
--k: the number of seed nodes to be selected.
--batch: the size of batch b selected each time. 
--seedfile: the file records the k seed nodes selected.
--time: the number of the algorithm repeated.

For example:
./exp_epic -dataset dataset/hep/ -model IC -epsilon 0.5 -k 500 -batch 50 -seedfile seed -time 1


