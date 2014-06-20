/*
 * BregmanFunctionalUnitTest.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "BregmanFunctionalUnitTest.hpp"

#include <Eigen/Dense>

#include "Minimizations/BregmanFunctional.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( BregmanFunctionalUnitTest );


void BregmanFunctionalUnitTest::setUp()
{
}


void BregmanFunctionalUnitTest::tearDown()
{
}

/**
 * Werte wurden wir folgt erzeugt:
 *
 * U=ones(2,2)-2.*rand(2,2)
 * alphas=ones(2,1)-2.*rand(2,1)
 * t=ones(2,10)-2.*rand(2,10)
 * X=ones(2,10)-2.*rand(2,10)
 * fval=zeros(1,10)
 * gval=zeros(2,10)
 * for i=1:10
 * 	[fval(:,i),gval(:,i)]=BregmanFunctional(t(:,i),X(:,i),U,alphas,inf,10.,1e-6)
 * endfor
 * fval
 * gval
 *
 * BregmanFunction ist aus Frank Schöpfer's Code und unverändert.
 */

void BregmanFunctionalUnitTest::oneoneNorm()
{
	const double p = 1.1;
	LpNorm lpnorm(p/(p-1.));
	LpNorm lpdualnorm(p);
	DualityMapping J_p(p);
	Eigen::MatrixXd U(2,2);
	U <<   -0.37133, -0.95779, 0.40013, 0.93263;
	Eigen::VectorXd alpha(2,1);
	alpha << 0.014897, -0.568410;
	Eigen::MatrixXd t(2,10);
	Eigen::MatrixXd X(2,10);
	t << -0.104574,0.559468,0.766960,0.148054,-0.923356,0.041415,0.015722,-0.033294,-0.570823,-0.051950,-0.914503,0.166395,-0.342347,-0.728231,-0.888305,0.778963,0.730667,-0.825560,-0.306347,-0.579328;
	X << 0.768586,-0.225475,0.470209,0.483057,-0.324979,-0.554573,0.030939,-0.352048,-0.103521,0.305334,-0.902129,-0.610977,0.396205,0.854957,0.562012,-0.658914,0.964865,-0.873672,0.717819,0.118345;
	{
		const double power = 1.01;
		Eigen::VectorXd expected(10);
		expected << 0.66508,0.99938,0.98106,1.99821,3.57436,1.10574,0.53693,1.67773,1.89644,1.21761;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << -0.057909,0.714264,-0.010482,-0.680474,-0.717641,0.717434,0.023131,-0.046084,-0.710026,-0.702523,-0.828367,1.136983,-0.540930,-2.262869,-2.362797,1.144847,-0.459392,-0.802523,-2.341643,-2.322302;
		BregmanFunctional d_p(lpdualnorm, J_p);
		for (size_t i=0; i<10; ++i) {
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t.col(i),X.col(i),U,alpha,power) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t.col(i),X.col(i),U,alpha,power) ) ) < 1e-4);
//			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << d_p.gradient(t.col(i),X.col(i),U,alpha,power).transpose() << ".\n";
			CPPUNIT_ASSERT( lpnorm(expected_gradient.col(i) - d_p.gradient(t.col(i),X.col(i),U,alpha,power) ) < 1e-4);
		}
	}
	{
		const double power = 1.1;
		Eigen::VectorXd expected(10);
		expected << 0.63198,0.91878,0.90229,1.93072,3.62365,1.03683,0.45599,1.59875,1.83592,1.13708;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << -0.0465225,0.7200278,-0.0099343,-0.7101301,-0.7956714,0.7459725,0.0231023,-0.0471788,-0.7470194,-0.6956763,-0.7877107,1.1510384,-0.5415234,-2.3351357,-2.5539356,1.2144422,-0.4597725,-0.8067273,-2.4321315,-2.3055647;
		BregmanFunctional d_p(lpdualnorm, J_p);
		for (size_t i=0; i<10; ++i) {
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t.col(i),X.col(i),U,alpha,power) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t.col(i),X.col(i),U,alpha,power) ) ) < 1e-4);
//			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << d_p.gradient(t.col(i),X.col(i),U,alpha,power).transpose() << ".\n";
			CPPUNIT_ASSERT( lpnorm(expected_gradient.col(i) - d_p.gradient(t.col(i),X.col(i),U,alpha,power) ) < 1e-4);
		}
	}
	{
		const double power = 2.;
		Eigen::VectorXd expected(10);
		expected << 0.529673,0.513800,0.513891,1.680977,5.231910,0.769151,0.047546,1.211317,1.676777,0.732566;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << 0.0036843,0.7803462,-0.0050663,-1.0859662,-2.2155091,1.1035507,0.0228205,-0.0592704,-1.2384149,-0.6307056,-0.6084441,1.2981238,-0.5467944,-3.2509634,-6.0318975,2.0864571,-0.4635044,-0.8531492,-3.6341332,-2.1467293;
		BregmanFunctional d_p(lpdualnorm, J_p);
		for (size_t i=0; i<10; ++i) {
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t.col(i),X.col(i),U,alpha,power) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t.col(i),X.col(i),U,alpha,power) ) ) < 1e-4);
//			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << d_p.gradient(t.col(i),X.col(i),U,alpha,power).transpose() << ".\n";
			CPPUNIT_ASSERT( lpnorm(expected_gradient.col(i) - d_p.gradient(t.col(i),X.col(i),U,alpha,power) ) < 1e-4);
		}
	}
}

void BregmanFunctionalUnitTest::twoNorm()
{
	const double p = 2.;
	LpNorm lpnorm(p/(p-1.));
	LpNorm lpdualnorm(p);
	DualityMapping J_p(p);
	Eigen::MatrixXd U(2,2);
	U <<   -0.37133, -0.95779, 0.40013, 0.93263;
	Eigen::VectorXd alpha(2,1);
	alpha << 0.014897, -0.568410;
	Eigen::MatrixXd t(2,10);
	Eigen::MatrixXd X(2,10);
	t << -0.104574,0.559468,0.766960,0.148054,-0.923356,0.041415,0.015722,-0.033294,-0.570823,-0.051950,-0.914503,0.166395,-0.342347,-0.728231,-0.888305,0.778963,0.730667,-0.825560,-0.306347,-0.579328;
	X << 0.768586,-0.225475,0.470209,0.483057,-0.324979,-0.554573,0.030939,-0.352048,-0.103521,0.305334,-0.902129,-0.610977,0.396205,0.854957,0.562012,-0.658914,0.964865,-0.873672,0.717819,0.118345;
	{
		const double power = 1.1;
		Eigen::VectorXd expected(10);
		expected << 0.62802,0.82295,0.71584,1.81901,2.79042,0.89188,0.28344,1.53771,1.45524,0.97226;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << -0.2744390,0.4635865,0.0070310,-0.4404360,-0.5786415,0.4809155,0.2166510,-0.3291053,-0.5252405,-0.4780703,-1.3188585,0.4904822,-0.5234721,-1.6400329,-2.0174474,0.5317536,-0.0138196,-1.4635940,-1.8695927,-1.7503391;
		BregmanFunctional d_p(lpdualnorm, J_p);
		for (size_t i=0; i<10; ++i) {
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t.col(i),X.col(i),U,alpha,power) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t.col(i),X.col(i),U,alpha,power) ) ) < 1e-4);
//			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << d_p.gradient(t.col(i),X.col(i),U,alpha,power).transpose() << ".\n";
			CPPUNIT_ASSERT( lpnorm(expected_gradient.col(i) - d_p.gradient(t.col(i),X.col(i),U,alpha,power) ) < 1e-4);
		}
	}
	{
		const double power = 2.;
		Eigen::VectorXd expected(10);
		expected << 0.52896,0.41386,0.38071,1.51650,3.19301,0.56201,-0.10539,1.14000,1.10983,0.59547;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << -0.0364139,0.4636278,0.0099964,-0.6344602,-1.2532194,0.6526892,0.1775279,-0.3778562,-0.7041339,-0.3567879,-0.7014935,0.4905796,-0.5404139,-2.0966682,-3.6643320,0.9372717,-0.1213626,-1.5904568,-2.3005447,-1.4595545;
		BregmanFunctional d_p(lpdualnorm, J_p);
		for (size_t i=0; i<10; ++i) {
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t.col(i),X.col(i),U,alpha,power) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t.col(i),X.col(i),U,alpha,power) ) ) < 1e-4);
//			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << d_p.gradient(t.col(i),X.col(i),U,alpha,power).transpose() << ".\n";
			CPPUNIT_ASSERT( lpnorm(expected_gradient.col(i) - d_p.gradient(t.col(i),X.col(i),U,alpha,power) ) < 1e-4);
		}
	}
	{
		const double power = 4.;
		Eigen::VectorXd expected(10);
		expected << 0.51837,0.16386,0.23653,1.62694,7.79111,0.56620,-0.31918,0.91932,1.05715,0.39978;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << 0.013798,0.463720,0.013185,-1.414160,-6.837598,1.295795,0.115627,-0.512366,-1.342928,-0.183546,-0.571259,0.490796,-0.558629,-3.931688,-17.297786,2.455497,-0.291517,-1.940487,-3.839391,-1.044191;
		BregmanFunctional d_p(lpdualnorm, J_p);
		for (size_t i=0; i<10; ++i) {
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t.col(i),X.col(i),U,alpha,power) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t.col(i),X.col(i),U,alpha,power) ) ) < 1e-4);
//			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << d_p.gradient(t.col(i),X.col(i),U,alpha,power).transpose() << ".\n";
			CPPUNIT_ASSERT( lpnorm(expected_gradient.col(i) - d_p.gradient(t.col(i),X.col(i),U,alpha,power) ) < 1e-4);
		}
	}
}

void BregmanFunctionalUnitTest::inftyNorm()
{
	// due to large/small values we look at relative precision
	const double p = LpNorm::Infinity;
	LpNorm lpnorm(p/(p-1.));
	LpNorm lpdualnorm(p);
	DualityMapping J_p(p);
	Eigen::MatrixXd U(2,2);
	U <<   -0.37133, -0.95779, 0.40013, 0.93263;
	Eigen::VectorXd alpha(2,1);
	alpha << 0.014897, -0.568410;
	Eigen::MatrixXd t(2,10);
	Eigen::MatrixXd X(2,10);
	t << -0.104574,0.559468,0.766960,0.148054,-0.923356,0.041415,0.015722,-0.033294,-0.570823,-0.051950,-0.914503,0.166395,-0.342347,-0.728231,-0.888305,0.778963,0.730667,-0.825560,-0.306347,-0.579328;
	X << 0.768586,-0.225475,0.470209,0.483057,-0.324979,-0.554573,0.030939,-0.352048,-0.103521,0.305334,-0.902129,-0.610977,0.396205,0.854957,0.562012,-0.658914,0.964865,-0.873672,0.717819,0.118345;
	{
		const double power = 1.1;
		Eigen::VectorXd expected(10);
		expected << 0.62787,0.81287,0.56263,1.81008,2.18416,0.87617,0.23439,1.53412,1.30917,0.92277;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << -0.29146,0.41462,0.35594,-0.40109,-0.40850,0.42878,0.37504,-0.36182,-0.39367,-0.37006,-1.35863,0.36329,0.31127,-1.53800,-1.55528,0.39627,0.36054,-1.54011,-1.52070,-1.46569;
		BregmanFunctional d_p(lpdualnorm, J_p);
		for (size_t i=0; i<10; ++i) {
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t.col(i),X.col(i),U,alpha,power) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t.col(i),X.col(i),U,alpha,power) )/expected(i) ) < 1e-4);
//			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << d_p.gradient(t.col(i),X.col(i),U,alpha,power).transpose() << ".\n";
			CPPUNIT_ASSERT( (lpnorm(expected_gradient.col(i) - d_p.gradient(t.col(i),X.col(i),U,alpha,power) ))/lpnorm(expected_gradient.col(i)) < 1e-4);
		}
	}
	{
		const double power = 2.;
		Eigen::VectorXd expected(10);
		expected << 0.52893,0.40382,0.29723,1.50379,2.03986,0.54061,-0.14379,1.13591,0.92445,0.55934;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << -0.039371,0.411033,0.173492,-0.575250,-0.689305,0.575867,0.288415,-0.414030,-0.478035,-0.256964,-0.708385,0.354918,-0.159334,-1.943942,-2.209787,0.739118,0.137097,-1.674771,-1.717351,-1.202072;
		BregmanFunctional d_p(lpdualnorm, J_p);
		for (size_t i=0; i<10; ++i) {
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t.col(i),X.col(i),U,alpha,power) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t.col(i),X.col(i),U,alpha,power) )/expected(i) ) < 1e-4);
//			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << d_p.gradient(t.col(i),X.col(i),U,alpha,power).transpose() << ".\n";
			CPPUNIT_ASSERT( (lpnorm(expected_gradient.col(i) - d_p.gradient(t.col(i),X.col(i),U,alpha,power) ))/lpnorm(expected_gradient.col(i)) < 1e-4);
		}
	}
	{
		const double power = 10.;
		Eigen::VectorXd expected(10);
		expected << 5.1825e-01,4.2115e-03,2.0604e-01,5.2869e+00,2.9000e+01,2.4914e+00,-4.1038e-01,8.9170e-01,9.7075e-01,3.3062e-01;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << 1.4897e-02,3.8049e-01,1.5072e-02,-1.3199e+01,-6.4800e+01,8.3874e+00,3.8602e-02,-1.3447e+00,-2.6001e+00,2.5507e-03,-5.6841e-01,2.8373e-01,-5.6796e-01,-3.1368e+01,-1.5164e+02,1.8947e+01,-5.0726e-01,-4.0753e+00,-6.6636e+00,-5.9719e-01;
		BregmanFunctional d_p(lpdualnorm, J_p);
		for (size_t i=0; i<10; ++i) {
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t.col(i),X.col(i),U,alpha,power) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t.col(i),X.col(i),U,alpha,power) )/expected(i) ) < 2e-4);
//			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << d_p.gradient(t.col(i),X.col(i),U,alpha,power).transpose() << ".\n";
			CPPUNIT_ASSERT( (lpnorm(expected_gradient.col(i) - d_p.gradient(t.col(i),X.col(i),U,alpha,power) ))/lpnorm(expected_gradient.col(i)) < 1e-4);
		}
	}
}
