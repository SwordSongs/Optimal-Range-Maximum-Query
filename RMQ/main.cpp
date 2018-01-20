/#include "rmq.h"
//#include <random>
//
//
//typedef util::rmq_fbst<int> rmq_fbst_type;
//typedef util::rmq_bender00<int> rmq_bender_type;
//typedef util::rmq_bender00_eco<int> rmq_bender_eco_type;
//typedef util::rmq_bender00_very_eco<int> rmq_bender_very_eco_type;
//typedef util::rmq_fischer11<int> rmq_fischer_type;
//
////---------------------------------------------------------------------------------------------------------------------
//
//struct inc_gen {
//	inc_gen() : curr(0) {}
//	int operator() () { return curr++; }
//	int curr;
//};
//
////---------------------------------------------------------------------------------------------------------------------
//
//int main()
//{
//	std::vector<int> input(10000);
//	std::generate(input.begin(), input.end(), inc_gen());
//	std::random_shuffle(input.begin(), input.end());
//
//	std::cout << "Constructing the FBST RMQ of size " << std::flush; begin_timing(); rmq_fbst_type rmq_fbst(input); rmq_fbst.init_sequential(); std::cout << rmq_fbst.size_in_bits() << " bits. " << std::flush; end_timing();
//	std::cout << "Constructing the standard variant Bender's RMQ of size " << std::flush; begin_timing(); rmq_bender_type rmq_bender(input); rmq_bender.init_sequential(); std::cout << rmq_bender.size_in_bits() << " bits. " << std::flush; end_timing();
//	std::cout << "Constructing the economic variant of Bender's RMQ of size " << std::flush; begin_timing(); rmq_bender_eco_type rmq_bender_eco(input); rmq_bender_eco.init_sequential(); std::cout << rmq_bender_eco.size_in_bits() << " bits. " << std::flush; end_timing();
//	std::cout << "Constructing the very economic variant of Bender's RMQ of size " << std::flush; begin_timing(); rmq_bender_very_eco_type rmq_bender_very_eco(input); rmq_bender_very_eco.init_sequential(); std::cout << rmq_bender_very_eco.size_in_bits() << " bits. " << std::flush; end_timing();
//	std::cout << "Constructing the Fischer's RMQ of size " << std::flush; begin_timing(); rmq_fischer_type rmq_fischer(input); rmq_fischer.init_sequential(); std::cout << rmq_fischer.size_in_bits() << " bits. " << std::flush; end_timing();
//
//	//-------------------------------------------------------------------------------------------------------------------------
//
//	//    std::cout << "Searching with FBST RMQ structure. " << std::flush;
//	//    begin_timing();
//	//    for(int i = 0; i < input.size(); ++i)
//	//    {
//	//        for(int j = i + 1; j < input.size(); ++j) { rmq_fbst.get_minimum(i, j); }
//	//        if((i + 1) % (input.size()/4) == 0) std::cout << " - " << std::flush;
//	//    }
//	//    end_timing();
//
//	//    //-------------------------------------------------------------------------------------------------------------------------
//
//	//    std::cout << "Searching with the standard variant of Bender's RMQ structure. " << std::flush;
//	//    begin_timing();
//	//    for(int i = 0; i < input.size(); ++i)
//	//    {
//	//        for(int j = i + 1; j < input.size(); ++j) { rmq_bender.get_minimum(i, j); }
//	//        if((i + 1) % (input.size()/4) == 0) std::cout << " - " << std::flush;
//	//    }
//	//    end_timing();
//
//	//    //-------------------------------------------------------------------------------------------------------------------------
//
//	//    std::cout << "Searching with the economic variant of Bender's RMQ structure. " << std::flush;
//	//    begin_timing();
//	//    for(int i = 0; i < input.size(); ++i)
//	//    {
//	//        for(int j = i + 1; j < input.size(); ++j) { rmq_bender_eco.get_minimum(i, j); }
//	//        if((i + 1) % (input.size()/4) == 0) std::cout << " - " << std::flush;
//	//    }
//	//    end_timing();
//
//	//    //-------------------------------------------------------------------------------------------------------------------------
//
//	//    std::cout << "Searching with the very economic variant of Bender's RMQ structure. " << std::flush;
//	//    begin_timing();
//	//    for(int i = 0; i < input.size(); ++i)
//	//    {
//	//        for(int j = i + 1; j < input.size(); ++j) { rmq_bender_very_eco.get_minimum(i, j); }
//	//        if((i + 1) % (input.size()/4) == 0) std::cout << " - " << std::flush;
//	//    }
//	//    end_timing();
//
//	//-------------------------------------------------------------------------------------------------------------------------
//
//	std::cout << "Searching with the Fischer's RMQ structure. " << std::flush;
//	begin_timing();
//	for (int i = 0; i < input.size(); ++i)
//	{
//		for (int j = i + 1; j < input.size(); ++j) { rmq_fischer.get_minimum(i, j); }
//		if ((i + 1) % (input.size() / 4) == 0) std::cout << " - " << std::flush;
//	}
//	end_timing();
//
//	//-------------------------------------------------------------------------------------------------------------------------
//
//	//    std::cout << "Searching via scanning. " << std::flush;
//	//    begin_timing();
//	//    for(int i = 0; i < input.size(); ++i)
//	//    {
//	//        for(int j = i + 1; j < input.size(); ++j) { util::get_minimum(i, j, input); }
//	//        if((i + 1) % (input.size()/4) == 0) std::cout << " - " << std::flush;
//	//    }
//	//    end_timing();
//
//	//========================================================================================================================
//
//	//    for(int i = 0; i < rmq_fbst.size(); ++i)
//	//    {
//	//        for(int j = i + 1; j < rmq_fbst.size(); ++j)
//	//        {
//	//            if(util::get_minimum(i, j, input) != rmq_fbst[rmq_fbst.get_minimum(i, j)])
//	//            {
//	//                std::cout << " error on " << i << "," << j << " => was "<< rmq_fbst.get_minimum(i, j) << " should be " << util::get_minimum(i, j, input) << std::endl;
//	//            } //else std::cout << "ok on " << i << "," << j << " =>" << std::endl;
//
//	//            if(i == rmq_fischer.size() - 2 && j == rmq_fischer.size() - 1) std::cout << "OK: Google" << std::endl;
//	//        }
//	//    }
//
//	//    for(int i = 0; i < rmq_fischer.size(); ++i)
//	//    {
//	//        for(int j = i + 1; j < rmq_fischer.size(); ++j)
//	//        {
//	//            if(util::get_minimum(i, j, input) != rmq_fischer[rmq_fischer.get_minimum(i, j)])
//	//            {
//	//                std::cout << " error on " << i << "," << j << " => was "<< rmq_fischer.get_minimum(i, j) << " should be " << util::get_minimum(i, j, input) << std::endl;
//	//            } //else std::cout << "ok on " << i << "," << j << " =>" << std::endl;
//
//	//            if(i == rmq_fischer.size() - 2 && j == rmq_fischer.size() - 1) std::cout << "OK: Fischer" << std::endl;
//	//        }
//	//    }
//
//	//    for(int i = 0; i < rmq_bender_very_eco.size(); ++i)
//	//    {
//	//        for(int j = i + 1; j < rmq_bender_very_eco.size(); ++j)
//	//        {
//	//            if(util::get_minimum(i, j, input) != rmq_bender_very_eco[rmq_bender_very_eco.get_minimum(i, j)])
//	//            {
//	//                std::cout << " error on " << i << "," << j << " => was "<< rmq_bender_very_eco.get_minimum(i, j) << " should be " << util::get_minimum(i, j, input) << std::endl;
//	//            }// else std::cout << "ok on " << i << "," << j << " =>" << std::endl;
//
//	//            if(i == rmq_fischer.size() - 2 && j == rmq_fischer.size() - 1) std::cout << "OK: Eco+ Bender" << std::endl;
//	//        }
//	//    }
//
//	//    for(int i = 0; i < rmq_bender_eco.size(); ++i)
//	//    {
//	//        for(int j = i + 1; j < rmq_bender_eco.size(); ++j)
//	//        {
//	//            if(util::get_minimum(i, j, input) != rmq_bender_eco[rmq_bender_eco.get_minimum(i, j)])
//	//            {
//	//                std::cout << " error on " << i << "," << j << " => was "<< rmq_bender_eco.get_minimum(i, j) << " should be " << util::get_minimum(i, j, input) << std::endl;
//	//            }// else std::cout << "ok on " << i << "," << j << " =>" << std::endl;
//
//	//            if(i == rmq_fischer.size() - 2 && j == rmq_fischer.size() - 1) std::cout << "OK: Eco Bender" << std::endl;
//	//        }
//	//    }
//
//	//    for(int i = 0; i < rmq_bender.size(); ++i)
//	//    {
//	//        for(int j = i + 1; j < rmq_bender.size(); ++j)
//	//        {
//	//            if(util::get_minimum(i, j, input) != rmq_bender[rmq_bender.get_minimum(i, j)])
//	//            {
//	//                std::cout << " error on " << i << "," << j << " => was "<< rmq_bender.get_minimum(i, j) << " should be " << util::get_minimum(i, j, input) << std::endl;
//	//            }// else std::cout << "ok on " << i << "," << j << " =>" << std::endl;
//
//	//            if(i == rmq_fischer.size() - 2 && j == rmq_fischer.size() - 1) std::cout << "OK: Bender" << std::endl;
//	//        }
//	//    }
//
//	return 0;
//}
