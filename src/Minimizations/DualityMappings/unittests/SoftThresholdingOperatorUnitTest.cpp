/*
 * SoftThresholdingOperatorUnitTest.cpp
 *
 *  Created on: Oct 13, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "SoftThresholdingOperatorUnitTest.hpp"

#include <Eigen/Dense>

#include "Minimizations/DualityMappings/SoftThresholdingOperator.hpp"
#include "Minimizations/MinimizationExceptions.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( SoftThresholdingOperatorUnitTest );


void SoftThresholdingOperatorUnitTest::setUp()
{
}


void SoftThresholdingOperatorUnitTest::tearDown()
{
}

/** We generate test vectors as follows via octave:
 *
 * x=ones(2,10)-2.*rand(2,10)
 * gval=zeros(2,10)
 * for i=1:10
 * 	gval(:,i)=SoftThresholdingOperator(x, p, power, 1e-6)
 * endfor
 * gval
 *
 *
 */

void SoftThresholdingOperatorUnitTest::oneNorm()
{
	SoftThresholdingOperator S;
	Eigen::VectorXd x(10);
	x << 0.204691,-0.799513,0.056042,0.364664,0.039179,-0.272607,-0.851628,0.720586,-0.058074,-0.529929;
	{
		const double lambda = .9;
		S.setLambda(lambda);
		Eigen::VectorXd expected(10);
		expected << 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.;
		const Eigen::VectorXd compare = S(x);
//			std::cout << "Expecting " << expected.transpose()
//					<< " and got " << compare.transpose() << ".\n";
		CPPUNIT_ASSERT( expected.isApprox(compare, 1e-4)  );
	}
	{
		const double lambda = .4;
		S.setLambda(lambda);
		Eigen::VectorXd expected(10);
		expected << 0.,-0.399513,0.,0.,0.,0.,-0.451628,0.320586,0.,-0.129929;
		const Eigen::VectorXd compare = S(x);
//			std::cout << "Expecting " << expected.transpose()
//					<< " and got " << compare.transpose() << ".\n";
		CPPUNIT_ASSERT( expected.isApprox(compare, 1e-4)  );
	}
	{
		const double lambda = .1;
		S.setLambda(lambda);
		Eigen::VectorXd expected(10);
		expected << 0.104691,-0.699513,0.,0.264664,0.,-0.172607,-0.751628,0.620586,0.,-0.429929;
		const Eigen::VectorXd compare = S(x);
//			std::cout << "Expecting " << expected.transpose()
//					<< " and got " << compare.transpose() << ".\n";
		CPPUNIT_ASSERT( expected.isApprox(compare, 1e-4)  );
	}
	{
		const double lambda = .001;
		S.setLambda(lambda);
		Eigen::VectorXd expected(10);
		expected << 0.203691,-0.798513,0.055042,0.363664,0.038179,-0.271607,-0.850628,0.719586,-0.057074,-0.528929;
		const Eigen::VectorXd compare = S(x);
//			std::cout << "Expecting " << expected.transpose()
//					<< " and got " << compare.transpose() << ".\n";
		CPPUNIT_ASSERT( expected.isApprox(compare, 1e-4)  );
	}
	{
		const double lambda = .0;
		S.setLambda(lambda);
		Eigen::VectorXd expected(x);
		const Eigen::VectorXd compare = S(x);
//			std::cout << "Expecting " << expected.transpose()
//					<< " and got " << compare.transpose() << ".\n";
		CPPUNIT_ASSERT( expected.isApprox(compare, 1e-4)  );
	}
}
