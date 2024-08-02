#include <iostream>
#include <string>
#include <ilcplex/ilocplex.h>

#include "src/timer.h"
#include "src/heuristic_solution.h"
#include "src/kernel_pump/kernel_pump.h"
#include "src/problem.h"

int main()
{
	try
	{
		// std::string problem_path = "/opt/ibm/ILOG/CPLEX_Studio2211/cplex/examples/data/p0033.mps";
		// std::string problem_path = "/home/lucas/Downloads/p0033.mps";
		// std::string problem_path = "/home/lucas/Downloads/mining.mps"; // KP achar melhor que cplex, mas demora. FP falha? (demora muito cada iteracao..nao deixei rodar tudo).
		// std::string problem_path = "/home/lucas/Downloads/hanoi5.mps"; // nao acha
		// std::string problem_path = "/home/lucas/Downloads/irp.mps";
		// std::string problem_path = "/home/lucas/Downloads/chromaticindex32-8.mps"; // acha igual para todos.
		// std::string problem_path = "/home/lucas/Documentos/Research/kernel-pump/codes/new_pb.lp";
		// std::string problem_path = "/home/lucas/Downloads/ramos3.mps"; // kp melhor que cplex e FP
		// std::string problem_path = "/home/lucas/Downloads/neos-578379.mps"; // nao acha (mas tem sol trivial zero!)
		// std::string problem_path = "/home/lucas/Downloads/reblock115.mps"; // feas pump foi melhor, mas KP melhor que cplex.
		// std::string problem_path = "/home/lucas/Downloads/supportcase30.mps"; // não acha (nenhum metodo)
		// std::string problem_path = "/home/lucas/Downloads/r-stop-dp-r50_4.lp"; // nao acha!
		// std::string problem_path = "/home/lucas/Downloads/r-stop-dp-c25_2.lp"; // FP e KP acharam melhor que o cplex. (mas so quando coloquei pra forcar adicao de bucket ao kernel enquanto ainda nao acha viavel)
		// std::string problem_path = "/home/lucas/Downloads/toll-like.mps"; // KP melhor que cplex... FP falhou.
		// std::string problem_path = "/home/lucas/Downloads/disctom.mps"; // achou a mesma para todos.
		// std::string problem_path = "/home/lucas/Downloads/graph40-20-1rand.mps"; // achou melhor que cplex (mas demorou mais). FP falhou.
		// std::string problem_path = "/home/lucas/Downloads/z26.mps"; // só CPLEX acha mas tem solucao trivial.
		// std::string problem_path = "/home/lucas/Downloads/stein15inf.mps"; // infeasible
		std::string problem_path = "/home/lucas/Downloads/rmine13.mps"; // achou melhor que o CPLEX (solucao trivial do cplex). FP falhou (demora muito).

		Timestamp *ti = NewTimestamp();
		Timer *timer = GetTimer();
		bool multithreading = false;
		bool reset_fp_initial_basis_at_new_loop = false;
		bool sort_by_fractional_part = true;
		bool always_force_bucket_vars_into_kernel = true;
		bool solve_kernel_pump = false;
		bool solve_feas_pump = false;
		bool solve_cplex = true;

		Problem pb(problem_path, multithreading);

		if (solve_kernel_pump)
		{
			std::cout << " *** KERNEL PUMP" << std::endl;
			pb.Reset();
			timer->Clock(ti);
			KernelPump kp;
			kp.Init(&pb);
			kp.Run(K_KS_MAX_SIZE_BUCKET, K_KS_MIN_TIME_LIMIT, K_KS_MAX_TIME_LIMIT, K_KS_DECAY_FACTOR_TIME_LIMIT, sort_by_fractional_part, reset_fp_initial_basis_at_new_loop, always_force_bucket_vars_into_kernel);
			std::cout << "time spent: " << timer->CurrentElapsedTime(ti) << " s" << std::endl;
		}

		if (solve_feas_pump)
		{
			std::cout << " *** FEAS PUMP" << std::endl;
			pb.Reset();
			timer->Clock(ti);
			FeasibilityPump fp;
			fp.Init(&pb);
			fp.Run(reset_fp_initial_basis_at_new_loop, 900);
			std::cout << "time spent: " << timer->CurrentElapsedTime(ti) << " s" << std::endl;
		}

		if (solve_cplex)
		{
			std::cout << " *** CPLEX" << std::endl;
			pb.Reset();
			timer->Clock(ti);
			pb.Solve(false, true);
			auto value_opt = pb.curr_obj_value();
			if (value_opt)
				std::cout << "first int solution value: " << *value_opt << std::endl;
			else
				std::cout << "no int solution" << std::endl;
			std::cout << "time spent: " << timer->CurrentElapsedTime(ti) << " s" << std::endl;
		}

		auto trivial_sol = boost::dynamic_bitset<>(pb.num_vars(), 0);
		pb.ComputeSolutionValue(trivial_sol);
		auto value_opt = pb.curr_obj_value();
		if (value_opt)
			std::cout << "trivial solution value: " << *value_opt << std::endl;
		else
			std::cout << "no int trivial solution" << std::endl;

		delete ti;
		ti = nullptr;
		DeleteTimer();
	}
	catch (const std::runtime_error &re)
	{
		std::cout << "Runtime error: " << re.what() << std::endl;
	}
	catch (const std::exception &ex)
	{
		std::cout << "Error occurred: " << ex.what() << std::endl;
	}
	catch (const int &error)
	{
		std::cout << "Error occurred: " << error << std::endl;
	}
	catch (IloException &e)
	{
		std::cout << "Concert Exception: " << e << std::endl;
	}
	catch (const char *e)
	{
		std::cout << e << std::endl;
	}
	catch (const std::string &e)
	{
		std::cout << e << std::endl;
	}
	catch (...)
	{
		std::cout << "Unknown failure occurred. Possible memory corruption" << std::endl;
	}

	return 0;
}
