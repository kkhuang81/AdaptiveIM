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
			
			const double i_max = ceil(log((2. + 2. * epsilon_a / 3)*left_n / epsilon_a / epsilon_a) / log(2)) + 1;
						
			const double ai = log(2 * i_max / delta);

			double sample = (log(2 / delta) + Math::logcnk(left_n, batch)) / batch;

			//debug
			const double delta_debug = 0.01*arg.epsilon/left_n;
			const double epsilon_debug = (arg.epsilon - delta_debug) / (1 - delta_debug);
			const double i_max_debug = ceil(log((2. + 2. * epsilon_debug / 3)*left_n / epsilon_debug / epsilon_debug) / log(2)) + 1;
			const double ai_debug = log(3 * i_max_debug / delta_debug);
			//debug

			cout << "sample: " << sample << "\ti_max:" << i_max << endl;
			vector<int>batch_set;

			for (int i = 0; i < i_max; i++)			
			{
				g.build_hyper_graph_r((size_t)sample);
				//cout << "size: " << (int64)(sample) << endl;
				batch_set.clear();
				//calculate the upper bound here
				double influence = g.build_seedset(batch, batch_set);  //we need to be careful here.							
							    
				double cover = g.Coverage(sample, batch_set);				

				double lower = sqr(sqrt(cover + 2. * ai / 9.) - sqrt(ai / 2.)) - ai / 18;				
				
				//debug
				double upper = sqr(sqrt(influence + ai / 2.) + sqrt(ai / 2.));				
				double lower_debug = sqr(sqrt(cover + 2. * ai_debug / 9.) - sqrt(ai_debug / 2.)) - ai_debug / 18;
				double upper_debug = sqr(sqrt(influence + ai_debug / 2.) + sqrt(ai_debug / 2.));
				//debug
				
				//if (g.hyperGT.size() != g.hyperGT_2.size())cout << "ALERT" << endl;
				double ratio = lower / influence;
				//cout << "cover " << cover << " lower " << lower << " influence " << influence << " ratio " << ratio << " sample " << sample << endl;
				//cout << influence << "\t" << g.hyperGT.size() << "\t" << lower << "\t" << g.hyperGT_2.size() << "\t" << ratio << endl;
				cout << "round: " << i << "\tsample: " << sample << endl;
				cout << "EPIC: lower/cover: " << lower / cover << "\tinfluence/upper: " << influence / upper << "\tcover/influence: " << cover / influence << "\tratio: " << ratio << "\tratio1: " << lower / upper << "\trequired_ratio: " << alpha*(1. - epsilon_prime) << endl;
				cout << "OPIM-C: lower/cover: " << lower_debug / cover << "\tinfluence/upper: " << influence / upper_debug << "\tcover/influence: " << cover / influence << "\tratio: " << lower_debug / influence << "\tratio1: " << lower_debug / upper_debug << "\trequired_ratio: " << alpha*(1. - epsilon_debug) << endl << endl;
				//cout << "round: " << i << "\tratio: " << ratio << "\tratio1: " << lower / upper << "\trequired: " << alpha*(1. - epsilon_prime) << endl << endl;

				if (ratio > alpha*(1.- epsilon_prime))
				{					
					//cout << "cover: " << cover << "\tlower: " << lower << "\tinfluence: " << influence << "\tupper: " << upper << "\tsample: " << sample << endl;
					//cout << "round: " << i << "\tratio: " << ratio << "\tratio1: " << lower / upper << "\trequired: " << alpha*(1. - epsilon_prime)<< endl << endl;
					rr_num += (g.hyperGT.size() + g.hyperGT_2.size());
					for (auto it : batch_set)g.seedSet.push_back(it);
					g.realization(batch_set);
					//exit(1);
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
			//if (arg.batch == 1)alpha = 1;

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
					const double delta = 0.01*arg.epsilon*arg.batch / left_n;
					//const double epsilon_prime = (arg.epsilon - factor*delta) / (1 - delta);
					const double epsilon_prime = 1.*(arg.batch*arg.epsilon - delta*left_n) / (arg.batch - delta*left_n);
					const double epsilon_a = epsilon_prime / (1 - epsilon_prime);
					//cout << "epsilon_prime: " << epsilon_prime << " epsilon_prime: " << epsilon_prime << endl;
					AdaptiveSelect(g, arg, alpha, epsilon_prime, delta, epsilon_a);
				}

				high_resolution_clock::time_point endTime = high_resolution_clock::now();
				duration<double> interval = duration_cast<duration<double>>(endTime - startTime);

				total_time += (double)interval.count();
				total_spread += (g.n - g.NumcurNode);

				cout << "SingleSpread " << g.n - g.NumcurNode << endl;
				cout << "SingleRuntime " << (double)interval.count() << endl;				
            }            
			//cout << "SampleTime(s) " << g.rand_time / arg.time << endl;
			cout << "RunningTime(s) " << total_time / arg.time << endl;
            disp_mem_usage();
			cout << "TotalSample " << rr_num / arg.time << endl;
			cout << "AdaptiveSpread " << total_spread / arg.time << endl;			
        }
};

unsigned long long  AIM::rr_num = 0;
unsigned int  AIM::left_n = 0;