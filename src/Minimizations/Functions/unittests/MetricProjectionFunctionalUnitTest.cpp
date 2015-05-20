/*
 * MetricProjectionFunctionalUnitTest.cpp
 *
 *  Created on: May 20, 2015
 *      Author: heber
 */

#include "MetricProjectionFunctionalUnitTest.hpp"

#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <iterator>
#include <limits>

#include "Log/Logging.hpp"

#include "Math/Helpers.hpp"
#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Functions/MetricProjectionFunctional.hpp"
#include "Minimizations/Mappings/Specifics/L1DualityMapping.hpp"
#include "Minimizations/Mappings/Specifics/LpDualityMapping.hpp"
#include "Minimizations/Mappings/Specifics/LInfinityDualityMapping.hpp"
#include "Minimizations/Norms/L1Norm.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Norms/LInfinityNorm.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( MetricProjectionFunctionalUnitTest );

using namespace boost::assign;

// static instances
const double MetricProjectionFunctionalUnitTest::tolerance = 1e-4;

void MetricProjectionFunctionalUnitTest::setUp()
{
}


void MetricProjectionFunctionalUnitTest::tearDown()
{
}

static std::vector<SpaceElement_ptr_t> getVectorOfSpaceElements(
		const NormedSpace_ptr_t &_dual_space,
		const Eigen::MatrixXd &_values)
{
	std::vector<SpaceElement_ptr_t> U;
	for (int i=0;i<_values.outerSize();++i)
		U.push_back(
				ElementCreator::create(_dual_space,_values.col(i)));
	assert( U.size() == (size_t)_values.outerSize() );
	return U;
}

static std::vector< std::vector<double> > getVectorOfVectors(
		const Eigen::MatrixXd &_values)
{
	std::vector< std::vector<double> > U(_values.outerSize(), std::vector<double>());
	for(int i=0;i<_values.outerSize(); ++i)
		for(int j=0;j<_values.innerSize(); ++j)
			U[i].push_back( _values.col(i)[j] );
	return U;
}

//void MetricProjectionFunctionalUnitTest::oneNorm()
//{
//	// due to large/small values we look at relative precision
//	const double p = Helpers::ConjugateValue(1.);
//	Eigen::MatrixXd Utemp(2,2);
//	Utemp << 0.68025,0.66385, -0.84814,-0.76287;
//	Eigen::MatrixXd ttemp(2,10);
//	ttemp << -0.091707,-0.897343,-0.127139,0.300883,-0.454444,0.967974,-0.262831,-0.668766,-0.670414,-0.455131,0.543059,0.865680,-0.300487,0.028029,-0.441100,0.019350,-0.775496,-0.908793,0.977006,-0.275997;
//	std::vector< std::vector<double> > t =
//			getVectorOfVectors(ttemp);
//	Eigen::MatrixXd Xtemp(2,10);
//	Xtemp << -0.379750,-0.943336,-0.883381,0.807463,-0.697762,-0.066113,0.894920,-0.081799,0.728086,0.113930,-0.690921,0.420590,0.819499,-0.726720,0.284505,0.338650,-0.116910,0.517283,0.074392,0.520657;
//	{
//		const double power = Helpers::ConjugateValue(1.1);
//		const NormedSpace_ptr_t SpaceX =
//				NormedSpaceFactory::createLpInstance(
//						Xtemp.innerSize(), p, power);
//		std::vector<SpaceElement_ptr_t> U =
//				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
//		const Mapping &J_1 = *SpaceX->getDualSpace()->getDualityMapping();
//		const MetricProjectionFunctional d_p(
//				*SpaceX->getDualSpace()->getNorm(),
//				dynamic_cast<const PowerTypeDualityMapping &>(J_1),
//				J_1.getPower(),U);
//		Eigen::VectorXd expected(10);
//		expected << 0.82176,1.01381,1.07100,0.91063,0.59857,1.76535,2.71796,1.92243,0.52698,0.70196;
//		Eigen::MatrixXd expected_gradient(2,10);
//		expected_gradient << -0.255682,1.472781,1.452913,-1.620815,-0.244916,1.543453,-1.763641,-1.700792,0.076645,-1.558400,-0.334443,1.221157,1.202610,-1.666661,-0.328093,1.287128,-1.799987,-1.741318,-0.138450,-1.608398;
//		for (size_t i=0; i<10; ++i) {
//			SpaceElement_ptr_t x =
//					ElementCreator::create(
//						SpaceX->getDualSpace(),
//						Xtemp.col(i));
////			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t[i],x) << ".\n";
//			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t[i],x) )/expected(i) ) < tolerance);
//			const std::vector<double> gval = d_p.gradient(t[i],x);
////			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << gval << ".\n";
//			Eigen::VectorXd difference_vector = expected_gradient.col(i);
//			for (int j=0;j<difference_vector.size();++j)
//				difference_vector[j] -= gval[j];
//			CPPUNIT_ASSERT( (difference_vector.norm())/(expected_gradient.col(i)).norm() < tolerance);
//		}
//	}
//	{
//		const double power = Helpers::ConjugateValue(2.);
//		const NormedSpace_ptr_t SpaceX =
//				NormedSpaceFactory::createLpInstance(
//						Xtemp.innerSize(), p, power);
//		std::vector<SpaceElement_ptr_t> U =
//				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
//		const Mapping &J_1 = *SpaceX->getDualSpace()->getDualityMapping();
//		const MetricProjectionFunctional d_p(
//				*SpaceX->getDualSpace()->getNorm(),
//				dynamic_cast<const PowerTypeDualityMapping &>(J_1),
//				J_1.getPower(),U);
//		Eigen::VectorXd expected(10);
//		expected << 0.41314,0.62818,0.66478,0.50207,0.28553,1.73849,3.38035,1.75062,0.13821,0.33762;
//		Eigen::MatrixXd expected_gradient(2,10);
//		expected_gradient << -0.260570,1.788871,1.563169,-1.668112,-0.176782,2.834724,-3.938726,-2.715546,0.044809,-1.130629,-0.337326,1.516221,1.305532,-1.710812,-0.287911,2.492506,-3.830390,-2.688573,-0.157226,-1.209081;
//		for (size_t i=0; i<10; ++i) {
//			SpaceElement_ptr_t x =
//					ElementCreator::create(
//						SpaceX->getDualSpace(),
//						Xtemp.col(i));
////			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t[i],x) << ".\n";
//			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t[i],x) )/expected(i) ) < tolerance);
//			const std::vector<double> gval = d_p.gradient(t[i],x);
////			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << gval << ".\n";
//			Eigen::VectorXd difference_vector = expected_gradient.col(i);
//			for (int j=0;j<difference_vector.size();++j)
//				difference_vector[j] -= gval[j];
//			CPPUNIT_ASSERT( (difference_vector.norm())/(expected_gradient.col(i)).norm() < tolerance);
//		}
//	}
//	{
//		const double power = Helpers::ConjugateValue(10.);
//		const NormedSpace_ptr_t SpaceX =
//				NormedSpaceFactory::createLpInstance(
//						Xtemp.innerSize(), p, power);
//		std::vector<SpaceElement_ptr_t> U =
//				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
//		const Mapping &J_1 = *SpaceX->getDualSpace()->getDualityMapping();
//		const MetricProjectionFunctional d_p(
//				*SpaceX->getDualSpace()->getNorm(),
//				dynamic_cast<const PowerTypeDualityMapping &>(J_1),
//				J_1.getPower(),U);
//		Eigen::VectorXd expected(10);
//		expected << 1.7741e-02,6.5155e-01,2.9734e-01,1.0730e-01,1.4355e-01,6.5141e+01,1.0328e+03,2.2887e+01,-1.6213e-01,1.0680e-01;
//		Eigen::MatrixXd expected_gradient(2,10);
//		expected_gradient << -3.1075e-01,9.5844e+00,2.9642e+00,-2.1581e+00,-8.7846e-02,5.2139e+02,-6.2629e+03,-2.0108e+02,-6.7892e-02,-1.3647e-01,-3.6692e-01,8.7932e+00,2.6134e+00,-2.1682e+00,-2.3546e-01,4.8655e+02,-5.8464e+03,-1.8786e+02,-2.2369e-01,-2.8105e-01;
//		for (size_t i=0; i<10; ++i) {
//			SpaceElement_ptr_t x =
//					ElementCreator::create(
//						SpaceX->getDualSpace(),
//						Xtemp.col(i));
////			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t[i],x) << ".\n";
//			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t[i],x) )/expected(i) ) < 4e-4);
//			const std::vector<double> gval = d_p.gradient(t[i],x);
////			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << gval << ".\n";
//			Eigen::VectorXd difference_vector = expected_gradient.col(i);
//			for (int j=0;j<difference_vector.size();++j)
//				difference_vector[j] -= gval[j];
//			CPPUNIT_ASSERT( (difference_vector.norm())/(expected_gradient.col(i)).norm() < tolerance);
//		}
//	}
//}
//
//void MetricProjectionFunctionalUnitTest::oneoneNorm()
//{
//	const double p = Helpers::ConjugateValue(1.1);
//	Eigen::MatrixXd Utemp(2,2);
//	Utemp <<   -0.37133, -0.95779, 0.40013, 0.93263;
//	Eigen::MatrixXd ttemp(2,10);
//	ttemp << -0.104574,0.559468,0.766960,0.148054,-0.923356,0.041415,0.015722,-0.033294,-0.570823,-0.051950,-0.914503,0.166395,-0.342347,-0.728231,-0.888305,0.778963,0.730667,-0.825560,-0.306347,-0.579328;
//	std::vector< std::vector<double> > t =
//			getVectorOfVectors(ttemp);
//	Eigen::MatrixXd Xtemp(2,10);
//	Xtemp << 0.768586,-0.225475,0.470209,0.483057,-0.324979,-0.554573,0.030939,-0.352048,-0.103521,0.305334,-0.902129,-0.610977,0.396205,0.854957,0.562012,-0.658914,0.964865,-0.873672,0.717819,0.118345;
//	{
//		const double power = Helpers::ConjugateValue(1.01);
//		const NormedSpace_ptr_t SpaceX =
//				NormedSpaceFactory::createLpInstance(
//						Xtemp.innerSize(), p, power);
//		std::vector<SpaceElement_ptr_t> U =
//				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
//		const Mapping &J_p = *SpaceX->getDualSpace()->getDualityMapping();
//		const MetricProjectionFunctional d_p(
//				*SpaceX->getDualSpace()->getNorm(),
//				dynamic_cast<const PowerTypeDualityMapping &>(J_p),
//				J_p.getPower(),U);
//		Eigen::VectorXd expected(10);
//		expected << 0.66508,0.99938,0.98106,1.99821,3.57436,1.10574,0.53693,1.67773,1.89644,1.21761;
//		Eigen::MatrixXd expected_gradient(2,10);
//		expected_gradient << -0.057909,0.714264,-0.010482,-0.680474,-0.717641,0.717434,0.023131,-0.046084,-0.710026,-0.702523,-0.828367,1.136983,-0.540930,-2.262869,-2.362797,1.144847,-0.459392,-0.802523,-2.341643,-2.322302;
//		for (size_t i=0; i<10; ++i) {
//			SpaceElement_ptr_t x =
//					ElementCreator::create(
//						SpaceX->getDualSpace(),
//						Xtemp.col(i));
////			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t[i],x) << ".\n";
//			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t[i],x) ) ) < tolerance);
//			const std::vector<double> gval = d_p.gradient(t[i],x);
////			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << gval << ".\n";
//			Eigen::VectorXd difference_vector = expected_gradient.col(i);
//			for (int j=0;j<difference_vector.size();++j)
//				difference_vector[j] -= gval[j];
//			CPPUNIT_ASSERT( difference_vector.norm() < tolerance);
//		}
//	}
//	{
//		const double power = Helpers::ConjugateValue(1.1);
//		const NormedSpace_ptr_t SpaceX =
//				NormedSpaceFactory::createLpInstance(
//						Xtemp.innerSize(), p, power);
//		std::vector<SpaceElement_ptr_t> U =
//				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
//		const Mapping &J_p = *SpaceX->getDualSpace()->getDualityMapping();
//		const MetricProjectionFunctional d_p(
//				*SpaceX->getDualSpace()->getNorm(),
//				dynamic_cast<const PowerTypeDualityMapping &>(J_p),
//				J_p.getPower(),U);
//		Eigen::VectorXd expected(10);
//		expected << 0.63198,0.91878,0.90229,1.93072,3.62365,1.03683,0.45599,1.59875,1.83592,1.13708;
//		Eigen::MatrixXd expected_gradient(2,10);
//		expected_gradient << -0.0465225,0.7200278,-0.0099343,-0.7101301,-0.7956714,0.7459725,0.0231023,-0.0471788,-0.7470194,-0.6956763,-0.7877107,1.1510384,-0.5415234,-2.3351357,-2.5539356,1.2144422,-0.4597725,-0.8067273,-2.4321315,-2.3055647;
//		for (size_t i=0; i<10; ++i) {
//			SpaceElement_ptr_t x =
//					ElementCreator::create(
//						SpaceX->getDualSpace(),
//						Xtemp.col(i));
////			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t[i],x) << ".\n";
//			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t[i],x) ) ) < tolerance);
//			const std::vector<double> gval = d_p.gradient(t[i],x);
////			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << gval << ".\n";
//			Eigen::VectorXd difference_vector = expected_gradient.col(i);
//			for (int j=0;j<difference_vector.size();++j)
//				difference_vector[j] -= gval[j];
//			CPPUNIT_ASSERT( difference_vector.norm() < tolerance);
//		}
//	}
//	{
//		const double power = Helpers::ConjugateValue(2.);
//		const NormedSpace_ptr_t SpaceX =
//				NormedSpaceFactory::createLpInstance(
//						Xtemp.innerSize(), p, power);
//		std::vector<SpaceElement_ptr_t> U =
//				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
//		const Mapping &J_p = *SpaceX->getDualSpace()->getDualityMapping();
//		const MetricProjectionFunctional d_p(
//				*SpaceX->getDualSpace()->getNorm(),
//				dynamic_cast<const PowerTypeDualityMapping &>(J_p),
//				J_p.getPower(),U);
//		Eigen::VectorXd expected(10);
//		expected << 0.529673,0.513800,0.513891,1.680977,5.231910,0.769151,0.047546,1.211317,1.676777,0.732566;
//		Eigen::MatrixXd expected_gradient(2,10);
//		expected_gradient << 0.0036843,0.7803462,-0.0050663,-1.0859662,-2.2155091,1.1035507,0.0228205,-0.0592704,-1.2384149,-0.6307056,-0.6084441,1.2981238,-0.5467944,-3.2509634,-6.0318975,2.0864571,-0.4635044,-0.8531492,-3.6341332,-2.1467293;
//		for (size_t i=0; i<10; ++i) {
//			SpaceElement_ptr_t x =
//					ElementCreator::create(
//						SpaceX->getDualSpace(),
//						Xtemp.col(i));
////			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t[i],x) << ".\n";
//			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t[i],x) ) ) < tolerance);
//			const std::vector<double> gval = d_p.gradient(t[i],x);
////			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << gval << ".\n";
//			Eigen::VectorXd difference_vector = expected_gradient.col(i);
//			for (int j=0;j<difference_vector.size();++j)
//				difference_vector[j] -= gval[j];
//			CPPUNIT_ASSERT( difference_vector.norm() < tolerance);
//		}
//	}
//}


/** Calculate reference values as follows:
 *
 * function:
 * for i=1:10
 * res(i)=(1/p)*norm(X(i,:)'-t(1,i)*U(:,1)-t(2,i)*U(:,2),p)^p
 * endfor
 *
 * gradient:
 *
 * for j=1:2
 * for i=1:10
 *  grad(j,i)=-1.*DualityMapping(X(i,:)'-t(1,i)*U(:,1)-t(2,i)*U(:,2),p,p,1e-8)'*U(:,j);
 * endfor
 * endfor
 */

void MetricProjectionFunctionalUnitTest::twoNorm()
{
	const double p = Helpers::ConjugateValue(2.); // changes nothing
	Eigen::MatrixXd Utemp(2,2);
	Utemp <<
-0.371330, -0.049584,
0.400130,  -0.046015;
	Eigen::MatrixXd ttemp(2,10);
	ttemp <<
0.864280,0.619754,0.898120,0.511959,0.965475,0.895577,0.739788,
0.175835,0.309747,0.750803,
0.111465,0.177248,0.912489,0.658087,0.345778,0.328143,0.098371,
0.899562,0.076553,0.047379;
	std::vector< std::vector<double> > t =
			getVectorOfVectors(ttemp);
	Eigen::MatrixXd Xtemp(2,10);
	Xtemp <<
0.696406,0.140523,0.685913,0.880953,0.990269,0.179758,0.282636,
0.704817,0.120739,0.590229,
0.288286,0.706717,0.323335,0.423850,0.346444,0.911419,0.010976,
0.276040,0.333290,0.933087;

	{
		const double power = Helpers::ConjugateValue(2.);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::createLpInstance(
						Xtemp.innerSize(), p, power);
		std::vector<SpaceElement_ptr_t> U =
				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
		const Mapping &J_p = *SpaceX->getDualSpace()->getDualityMapping();
		const MetricProjectionFunctional d_p(
				*SpaceX->getDualSpace()->getNorm(),
				dynamic_cast<const PowerTypeDualityMapping &>(J_p),
				J_p.getPower(),U);
		Eigen::VectorXd expected(10);
		expected <<
0.524501, 0.180983, 0.566765, 0.640136, 0.933161, 0.301109, 0.197388,
0.362403, 0.051350, 0.581162;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient <<
0.4007913,-0.0459179, 0.3929551, 0.3100879, 0.5167961,-0.0310635,
0.3210087, 0.2036646, 0.0037761, 0.0695452,
0.0483064, 0.0402986, 0.0530644, 0.0661964, 0.0666257, 0.0523540,
0.0149695, 0.0517663, 0.0216735, 0.0724192;
				for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x =
					ElementCreator::create(
						SpaceX->getDualSpace(),
						Xtemp.col(i));
//			std::cout << "Calculating (1/p)*norm([" << x << "] - "
//					<< t[i][0] << "*[" << U[0] << "]+" << t[i][1] << "*[" << U[1] << "], p)" << std::endl;
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t[i],x) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t[i],x) ) ) < tolerance);
			const std::vector<double> gval = d_p.gradient(t[i],x);
//			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << gval << ".\n";
			Eigen::VectorXd difference_vector = expected_gradient.col(i);
			for (int j=0;j<difference_vector.size();++j)
				difference_vector[j] -= gval[j];
			CPPUNIT_ASSERT( difference_vector.norm() < tolerance);
		}
	}
}

//void MetricProjectionFunctionalUnitTest::inftyNorm()
//{
//	// due to large/small values we look at relative precision
//	const double p = Helpers::ConjugateValue(std::numeric_limits<double>::infinity());
//	Eigen::MatrixXd Utemp(2,2);
//	Utemp <<   -0.37133, -0.95779, 0.40013, 0.93263;
//	Eigen::MatrixXd ttemp(2,10);
//	ttemp << -0.104574,0.559468,0.766960,0.148054,-0.923356,0.041415,0.015722,-0.033294,-0.570823,-0.051950,-0.914503,0.166395,-0.342347,-0.728231,-0.888305,0.778963,0.730667,-0.825560,-0.306347,-0.579328;
//	std::vector< std::vector<double> > t =
//			getVectorOfVectors(ttemp);
//	Eigen::MatrixXd Xtemp(2,10);
//	Xtemp << 0.768586,-0.225475,0.470209,0.483057,-0.324979,-0.554573,0.030939,-0.352048,-0.103521,0.305334,-0.902129,-0.610977,0.396205,0.854957,0.562012,-0.658914,0.964865,-0.873672,0.717819,0.118345;
//	{
//		const double power = Helpers::ConjugateValue(1.1);
//		const NormedSpace_ptr_t SpaceX =
//				NormedSpaceFactory::createLpInstance(
//						Xtemp.innerSize(), p, power);
//		std::vector<SpaceElement_ptr_t> U =
//				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
//		const Mapping &J_infty = *SpaceX->getDualSpace()->getDualityMapping();
//		const MetricProjectionFunctional d_p(
//				*SpaceX->getDualSpace()->getNorm(),
//				dynamic_cast<const PowerTypeDualityMapping &>(J_infty),
//				J_infty.getPower(),U);
//		Eigen::VectorXd expected(10);
//		expected << 0.62787,0.81287,0.56263,1.81008,2.18416,0.87617,0.23439,1.53412,1.30917,0.92277;
//		Eigen::MatrixXd expected_gradient(2,10);
//		expected_gradient << -0.29146,0.41462,0.35594,-0.40109,-0.40850,0.42878,0.37504,-0.36182,-0.39367,-0.37006,-1.35863,0.36329,0.31127,-1.53800,-1.55528,0.39627,0.36054,-1.54011,-1.52070,-1.46569;
//		for (size_t i=0; i<10; ++i) {
//			SpaceElement_ptr_t x =
//					ElementCreator::create(
//						SpaceX->getDualSpace(),
//						Xtemp.col(i));
////			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t[i],x) << ".\n";
//			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t[i],x) )/expected(i) ) < tolerance);
//			const std::vector<double> gval = d_p.gradient(t[i],x);
////			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << gval << ".\n";
//			Eigen::VectorXd difference_vector = expected_gradient.col(i);
//			for (int j=0;j<difference_vector.size();++j)
//				difference_vector[j] -= gval[j];
//			CPPUNIT_ASSERT( (difference_vector.norm())/(expected_gradient.col(i)).norm() < tolerance);
//		}
//	}
//	{
//		const double power = Helpers::ConjugateValue(2.);
//		const NormedSpace_ptr_t SpaceX =
//				NormedSpaceFactory::createLpInstance(
//						Xtemp.innerSize(), p, power);
//		std::vector<SpaceElement_ptr_t> U =
//				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
//		const Mapping &J_infty = *SpaceX->getDualSpace()->getDualityMapping();
//		const MetricProjectionFunctional d_p(
//				*SpaceX->getDualSpace()->getNorm(),
//				dynamic_cast<const PowerTypeDualityMapping &>(J_infty),
//				J_infty.getPower(),U);
//		Eigen::VectorXd expected(10);
//		expected << 0.52893,0.40382,0.29723,1.50379,2.03986,0.54061,-0.14379,1.13591,0.92445,0.55934;
//		Eigen::MatrixXd expected_gradient(2,10);
//		expected_gradient << -0.039371,0.411033,0.173492,-0.575250,-0.689305,0.575867,0.288415,-0.414030,-0.478035,-0.256964,-0.708385,0.354918,-0.159334,-1.943942,-2.209787,0.739118,0.137097,-1.674771,-1.717351,-1.202072;
//		for (size_t i=0; i<10; ++i) {
//			SpaceElement_ptr_t x =
//					ElementCreator::create(
//						SpaceX->getDualSpace(),
//						Xtemp.col(i));
////			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t[i],x) << ".\n";
//			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t[i],x) )/expected(i) ) < tolerance);
//			const std::vector<double> gval = d_p.gradient(t[i],x);
////			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << gval << ".\n";
//			Eigen::VectorXd difference_vector = expected_gradient.col(i);
//			for (int j=0;j<difference_vector.size();++j)
//				difference_vector[j] -= gval[j];
//			CPPUNIT_ASSERT( (difference_vector.norm())/(expected_gradient.col(i)).norm() < tolerance);
//		}
//	}
//	{
//		const double power = Helpers::ConjugateValue(10.);
//		const NormedSpace_ptr_t SpaceX =
//				NormedSpaceFactory::createLpInstance(
//						Xtemp.innerSize(), p, power);
//		std::vector<SpaceElement_ptr_t> U =
//				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
//		const Mapping &J_infty = *SpaceX->getDualSpace()->getDualityMapping();
//		const MetricProjectionFunctional d_p(
//				*SpaceX->getDualSpace()->getNorm(),
//				dynamic_cast<const PowerTypeDualityMapping &>(J_infty),
//				J_infty.getPower(),U);
//		Eigen::VectorXd expected(10);
//		expected << 5.1825e-01,4.2115e-03,2.0604e-01,5.2869e+00,2.9000e+01,2.4914e+00,-4.1038e-01,8.9170e-01,9.7075e-01,3.3062e-01;
//		Eigen::MatrixXd expected_gradient(2,10);
//		expected_gradient << 1.4897e-02,3.8049e-01,1.5072e-02,-1.3199e+01,-6.4800e+01,8.3874e+00,3.8602e-02,-1.3447e+00,-2.6001e+00,2.5507e-03,-5.6841e-01,2.8373e-01,-5.6796e-01,-3.1368e+01,-1.5164e+02,1.8947e+01,-5.0726e-01,-4.0753e+00,-6.6636e+00,-5.9719e-01;
//		for (size_t i=0; i<10; ++i) {
//			SpaceElement_ptr_t x =
//					ElementCreator::create(
//						SpaceX->getDualSpace(),
//						Xtemp.col(i));
////			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t[i],x) << ".\n";
//			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t[i],x) )/expected(i) ) < 2e-4);
//			const std::vector<double> gval = d_p.gradient(t[i],x);
////			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << gval << ".\n";
//			Eigen::VectorXd difference_vector = expected_gradient.col(i);
//			for (int j=0;j<difference_vector.size();++j)
//				difference_vector[j] -= gval[j];
//			CPPUNIT_ASSERT( (difference_vector.norm())/(expected_gradient.col(i)).norm() < tolerance);
//		}
//	}
//}
