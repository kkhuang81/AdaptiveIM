#include <chrono>
#include <ctime>
#include <ratio>

#define e exp(1)
#define c 2*(exp(1)-2)

using namespace std::chrono;

class Math{
    public:
        static double log2(int n){
            return log(n) / log(2);
        }
        static double logcnk(int n, int k) {
            double ans = 0;
            for (int i = n - k + 1; i <= n; i++)
            {
                ans += log(i);
            }
            for (int i = 1; i <= k; i++)
            {
                ans -= log(i);
            }
            return ans;
        }
};

class AIM
{
    private:       
        static unsigned long long rr_num;
		//static unsigned int active_node; //this is to record the number of the active node.
		static unsigned int left_n;		

		static void AdaptiveSelect(InfGraph &g, const Argument & arg, const double alpha, const double epsilon_prime, const double delta, const double epsilon_a)
		{			
			const unsigned int batch = arg.batch;			
			
			const double i_max = ceil(log((2. + 2 * epsilon_a / 3)*left_n / epsilon_a / epsilon_a) / log(2)) + 1;
						
			const double ai = log(2 * i_max / delta);

			double sample = (log(2 / delta) + Math::logcnk(left_n, batch)) / batch;
			
			vector<int>batch_set;

			for (int i = 0; i < i_max; i++)			
			{
				//Generate R1 and derive the seed set
				g.build_hyper_graph_r((int64)sample);
				
				batch_set.clear();
				//calculate the upper bound here
				double influence = g.build_seedset(batch, batch_set);  //we need to be careful here.							
							    
				double cover = g.Coverage(sample, batch_set);				

				double lower = sqr(sqrt(cover + 2. * ai / 9.) - sqrt(ai / 2.)) - ai / 18;				
								
				double ratio = lower / influence;				

				if (ratio > alpha*(1. - epsilon_prime))
				{				
					rr_num += (g.hyperGT.size() + g.hyperGT_2.size());
					for (auto it : batch_set)g.seedSet.push_back(it);
					g.realization(batch_set);				
					return;
				}				
				sample *= 2;
			}
			batch_set.clear();
			g.build_seedset(batch, batch_set);
			rr_num += (g.hyperGT.size() + g.hyperGT_2.size());
			for (auto it : batch_set)g.seedSet.push_back(it);
			g.realization(batch_set);			
		}
						
public:         	

	static void InfluenceMaximize(InfGraph &g, const Argument &arg)
        {                  			            								
			double total_spread = 0;			
			double total_time = 0;			
			
			double alpha = 1. - pow(1. - 1. / arg.batch, arg.batch);			

			const size_t rd = arg.k / arg.batch;
			//const double fdelta = 0.01;				//make the final delta const as 0.01;
			const double fdelta = 1. / g.n;
			const double epsilon_ast = arg.epsilon - sqrt(1. / (2 * rd)*log(1. / fdelta));
			cout << "epsilon_ast: " << epsilon_ast << endl;
			ASSERT(epsilon_ast > 0);			
			//cout << "factor: " << factor << endl;
			for (int i = 0; i < arg.time; i++)
            {				
				//the preparation work before each time.						
				g.load_possible_world(to_string(i), arg);
				g.init_hyper_graph();				
				//active_node = 0;
				left_n = g.n;

				high_resolution_clock::time_point startTime = high_resolution_clock::now();				

				while (g.seedSet.size() < arg.k)  //set the eta as 1;
 				{							
					left_n = g.NumcurNode;
					const double delta = 0.01*epsilon_ast*arg.batch / left_n;					
					const double epsilon_prime = 1.*(arg.batch*epsilon_ast - delta*left_n) / (arg.batch - delta*left_n);
					const double epsilon_a = epsilon_prime / (1 - epsilon_prime);					
					AdaptiveSelect(g, arg, alpha, epsilon_prime, delta, epsilon_a);
				}

				high_resolution_clock::time_point endTime = high_resolution_clock::now();
				duration<double> interval = duration_cast<duration<double>>(endTime - startTime);

				total_time += (double)interval.count();
				total_spread += (g.n - g.NumcurNode);

				cout << "SingleSpread " << g.n - g.NumcurNode << endl;
				cout << "SingleRuntime " << (double)interval.count() << endl;				
            }            			
			cout << "RunningTime(s) " << total_time / arg.time << endl;
            disp_mem_usage();
			cout << "TotalSample " << rr_num / arg.time << endl;
			cout << "AdaptiveSpread " << total_spread / arg.time << endl;			
        }
};

unsigned long long  AIM::rr_num = 0;
unsigned int  AIM::left_n = 0;
