/*
 * BregmanProjectionFunctionalUnitTest.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "BregmanProjectionFunctionalUnitTest.hpp"

#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <iterator>
#include <limits>

#include "Math/Helpers.hpp"
#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Functions/BregmanProjectionFunctional.hpp"
#include "Minimizations/Mappings/Specifics/L1DualityMapping.hpp"
#include "Minimizations/Mappings/Specifics/LpDualityMapping.hpp"
#include "Minimizations/Mappings/Specifics/LInfinityDualityMapping.hpp"
#include "Minimizations/Norms/L1Norm.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Norms/LInfinityNorm.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( BregmanProjectionFunctionalUnitTest );

using namespace boost::assign;

// static instances
const double BregmanProjectionFunctionalUnitTest::tolerance = 1e-4;

void BregmanProjectionFunctionalUnitTest::setUp()
{
}


void BregmanProjectionFunctionalUnitTest::tearDown()
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
 * BregmanFunctional ist aus Frank Schöpfer's Code und unverändert.
 */

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

template<class T>
std::ostream & operator<<(std::ostream &ost, const std::vector<T> & values)
{
	std::copy(
			values.begin(), values.end(),
			std::ostream_iterator<T>(ost, " "));
	return ost;
}

void BregmanProjectionFunctionalUnitTest::oneNorm()
{
	// due to large/small values we look at relative precision
	const double p = Helpers::ConjugateValue(1.);
	Eigen::MatrixXd Utemp(2,2);
	Utemp << 0.68025,0.66385, -0.84814,-0.76287;
	std::vector<double> alpha;
	alpha += -0.087261, -0.235115;
	Eigen::MatrixXd ttemp(2,10);
	ttemp << -0.091707,-0.897343,-0.127139,0.300883,-0.454444,0.967974,-0.262831,-0.668766,-0.670414,-0.455131,0.543059,0.865680,-0.300487,0.028029,-0.441100,0.019350,-0.775496,-0.908793,0.977006,-0.275997;
	std::vector< std::vector<double> > t =
			getVectorOfVectors(ttemp);
	Eigen::MatrixXd Xtemp(2,10);
	Xtemp << -0.379750,-0.943336,-0.883381,0.807463,-0.697762,-0.066113,0.894920,-0.081799,0.728086,0.113930,-0.690921,0.420590,0.819499,-0.726720,0.284505,0.338650,-0.116910,0.517283,0.074392,0.520657;
	{
		const double power = Helpers::ConjugateValue(1.1);
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						Xtemp.innerSize(), "lp", args);
		std::vector<SpaceElement_ptr_t> U =
				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
		const Mapping &J_1 = *SpaceX->getDualSpace()->getDualityMapping();
		const BregmanProjectionFunctional d_p(
				*SpaceX->getDualSpace()->getNorm(),
				dynamic_cast<const PowerTypeDualityMapping &>(J_1),
				J_1.getPower(),U,alpha);
		Eigen::VectorXd expected(10);
		expected << 0.82176,1.01381,1.07100,0.91063,0.59857,1.76535,2.71796,1.92243,0.52698,0.70196;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << -0.255682,1.472781,1.452913,-1.620815,-0.244916,1.543453,-1.763641,-1.700792,0.076645,-1.558400,-0.334443,1.221157,1.202610,-1.666661,-0.328093,1.287128,-1.799987,-1.741318,-0.138450,-1.608398;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x =
					ElementCreator::create(
						SpaceX->getDualSpace(),
						Xtemp.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t[i],x) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t[i],x) )/expected(i) ) < tolerance);
			const std::vector<double> gval = d_p.gradient(t[i],x);
//			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << gval << ".\n";
			Eigen::VectorXd difference_vector = expected_gradient.col(i);
			for (int j=0;j<difference_vector.size();++j)
				difference_vector[j] -= gval[j];
			CPPUNIT_ASSERT( (difference_vector.norm())/(expected_gradient.col(i)).norm() < tolerance);
		}
	}
	{
		const double power = Helpers::ConjugateValue(2.);
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						Xtemp.innerSize(), "lp", args);
		std::vector<SpaceElement_ptr_t> U =
				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
		const Mapping &J_1 = *SpaceX->getDualSpace()->getDualityMapping();
		const BregmanProjectionFunctional d_p(
				*SpaceX->getDualSpace()->getNorm(),
				dynamic_cast<const PowerTypeDualityMapping &>(J_1),
				J_1.getPower(),U,alpha);
		Eigen::VectorXd expected(10);
		expected << 0.41314,0.62818,0.66478,0.50207,0.28553,1.73849,3.38035,1.75062,0.13821,0.33762;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << -0.260570,1.788871,1.563169,-1.668112,-0.176782,2.834724,-3.938726,-2.715546,0.044809,-1.130629,-0.337326,1.516221,1.305532,-1.710812,-0.287911,2.492506,-3.830390,-2.688573,-0.157226,-1.209081;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x =
					ElementCreator::create(
						SpaceX->getDualSpace(),
						Xtemp.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t[i],x) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t[i],x) )/expected(i) ) < tolerance);
			const std::vector<double> gval = d_p.gradient(t[i],x);
//			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << gval << ".\n";
			Eigen::VectorXd difference_vector = expected_gradient.col(i);
			for (int j=0;j<difference_vector.size();++j)
				difference_vector[j] -= gval[j];
			CPPUNIT_ASSERT( (difference_vector.norm())/(expected_gradient.col(i)).norm() < tolerance);
		}
	}
	{
		const double power = Helpers::ConjugateValue(10.);
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						Xtemp.innerSize(), "lp", args);
		std::vector<SpaceElement_ptr_t> U =
				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
		const Mapping &J_1 = *SpaceX->getDualSpace()->getDualityMapping();
		const BregmanProjectionFunctional d_p(
				*SpaceX->getDualSpace()->getNorm(),
				dynamic_cast<const PowerTypeDualityMapping &>(J_1),
				J_1.getPower(),U,alpha);
		Eigen::VectorXd expected(10);
		expected << 1.7741e-02,6.5155e-01,2.9734e-01,1.0730e-01,1.4355e-01,6.5141e+01,1.0328e+03,2.2887e+01,-1.6213e-01,1.0680e-01;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << -3.1075e-01,9.5844e+00,2.9642e+00,-2.1581e+00,-8.7846e-02,5.2139e+02,-6.2629e+03,-2.0108e+02,-6.7892e-02,-1.3647e-01,-3.6692e-01,8.7932e+00,2.6134e+00,-2.1682e+00,-2.3546e-01,4.8655e+02,-5.8464e+03,-1.8786e+02,-2.2369e-01,-2.8105e-01;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x =
					ElementCreator::create(
						SpaceX->getDualSpace(),
						Xtemp.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t[i],x) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t[i],x) )/expected(i) ) < 4e-4);
			const std::vector<double> gval = d_p.gradient(t[i],x);
//			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << gval << ".\n";
			Eigen::VectorXd difference_vector = expected_gradient.col(i);
			for (int j=0;j<difference_vector.size();++j)
				difference_vector[j] -= gval[j];
			CPPUNIT_ASSERT( (difference_vector.norm())/(expected_gradient.col(i)).norm() < tolerance);
		}
	}
}

void BregmanProjectionFunctionalUnitTest::oneoneNorm()
{
	const double p = Helpers::ConjugateValue(1.1);
	Eigen::MatrixXd Utemp(2,2);
	Utemp <<   -0.37133, -0.95779, 0.40013, 0.93263;
	std::vector<double> alpha;
	alpha += 0.014897, -0.568410;
	Eigen::MatrixXd ttemp(2,10);
	ttemp << -0.104574,0.559468,0.766960,0.148054,-0.923356,0.041415,0.015722,-0.033294,-0.570823,-0.051950,-0.914503,0.166395,-0.342347,-0.728231,-0.888305,0.778963,0.730667,-0.825560,-0.306347,-0.579328;
	std::vector< std::vector<double> > t =
			getVectorOfVectors(ttemp);
	Eigen::MatrixXd Xtemp(2,10);
	Xtemp << 0.768586,-0.225475,0.470209,0.483057,-0.324979,-0.554573,0.030939,-0.352048,-0.103521,0.305334,-0.902129,-0.610977,0.396205,0.854957,0.562012,-0.658914,0.964865,-0.873672,0.717819,0.118345;
	{
		const double power = Helpers::ConjugateValue(1.01);
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						Xtemp.innerSize(), "lp", args);
		std::vector<SpaceElement_ptr_t> U =
				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
		const Mapping &J_p = *SpaceX->getDualSpace()->getDualityMapping();
		const BregmanProjectionFunctional d_p(
				*SpaceX->getDualSpace()->getNorm(),
				dynamic_cast<const PowerTypeDualityMapping &>(J_p),
				J_p.getPower(),U,alpha);
		Eigen::VectorXd expected(10);
		expected << 0.66508,0.99938,0.98106,1.99821,3.57436,1.10574,0.53693,1.67773,1.89644,1.21761;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << -0.057909,0.714264,-0.010482,-0.680474,-0.717641,0.717434,0.023131,-0.046084,-0.710026,-0.702523,-0.828367,1.136983,-0.540930,-2.262869,-2.362797,1.144847,-0.459392,-0.802523,-2.341643,-2.322302;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x =
					ElementCreator::create(
						SpaceX->getDualSpace(),
						Xtemp.col(i));
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
	{
		const double power = Helpers::ConjugateValue(1.1);
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						Xtemp.innerSize(), "lp", args);
		std::vector<SpaceElement_ptr_t> U =
				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
		const Mapping &J_p = *SpaceX->getDualSpace()->getDualityMapping();
		const BregmanProjectionFunctional d_p(
				*SpaceX->getDualSpace()->getNorm(),
				dynamic_cast<const PowerTypeDualityMapping &>(J_p),
				J_p.getPower(),U,alpha);
		Eigen::VectorXd expected(10);
		expected << 0.63198,0.91878,0.90229,1.93072,3.62365,1.03683,0.45599,1.59875,1.83592,1.13708;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << -0.0465225,0.7200278,-0.0099343,-0.7101301,-0.7956714,0.7459725,0.0231023,-0.0471788,-0.7470194,-0.6956763,-0.7877107,1.1510384,-0.5415234,-2.3351357,-2.5539356,1.2144422,-0.4597725,-0.8067273,-2.4321315,-2.3055647;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x =
					ElementCreator::create(
						SpaceX->getDualSpace(),
						Xtemp.col(i));
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
	{
		const double power = Helpers::ConjugateValue(2.);
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						Xtemp.innerSize(), "lp", args);
		std::vector<SpaceElement_ptr_t> U =
				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
		const Mapping &J_p = *SpaceX->getDualSpace()->getDualityMapping();
		const BregmanProjectionFunctional d_p(
				*SpaceX->getDualSpace()->getNorm(),
				dynamic_cast<const PowerTypeDualityMapping &>(J_p),
				J_p.getPower(),U,alpha);
		Eigen::VectorXd expected(10);
		expected << 0.529673,0.513800,0.513891,1.680977,5.231910,0.769151,0.047546,1.211317,1.676777,0.732566;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << 0.0036843,0.7803462,-0.0050663,-1.0859662,-2.2155091,1.1035507,0.0228205,-0.0592704,-1.2384149,-0.6307056,-0.6084441,1.2981238,-0.5467944,-3.2509634,-6.0318975,2.0864571,-0.4635044,-0.8531492,-3.6341332,-2.1467293;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x =
					ElementCreator::create(
						SpaceX->getDualSpace(),
						Xtemp.col(i));
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

void BregmanProjectionFunctionalUnitTest::twoNorm()
{
	const double p = Helpers::ConjugateValue(2.); // changes nothing
	Eigen::MatrixXd Utemp(2,2);
	Utemp <<   -0.37133, -0.95779, 0.40013, 0.93263;
	std::vector<double> alpha;
	alpha += 0.014897, -0.568410;
	Eigen::MatrixXd ttemp(2,10);
	ttemp << -0.104574,0.559468,0.766960,0.148054,-0.923356,0.041415,0.015722,-0.033294,-0.570823,-0.051950,-0.914503,0.166395,-0.342347,-0.728231,-0.888305,0.778963,0.730667,-0.825560,-0.306347,-0.579328;
	std::vector< std::vector<double> > t =
			getVectorOfVectors(ttemp);
	Eigen::MatrixXd Xtemp(2,10);
	Xtemp << 0.768586,-0.225475,0.470209,0.483057,-0.324979,-0.554573,0.030939,-0.352048,-0.103521,0.305334,-0.902129,-0.610977,0.396205,0.854957,0.562012,-0.658914,0.964865,-0.873672,0.717819,0.118345;
	{
		const double power = Helpers::ConjugateValue(1.1);
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						Xtemp.innerSize(), "lp", args);
		std::vector<SpaceElement_ptr_t> U =
				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
		const Mapping &J_p = *SpaceX->getDualSpace()->getDualityMapping();
		const BregmanProjectionFunctional d_p(
				*SpaceX->getDualSpace()->getNorm(),
				dynamic_cast<const PowerTypeDualityMapping &>(J_p),
				J_p.getPower(),U,alpha);
		Eigen::VectorXd expected(10);
		expected << 0.62802,0.82295,0.71584,1.81901,2.79042,0.89188,0.28344,1.53771,1.45524,0.97226;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << -0.2744390,0.4635865,0.0070310,-0.4404360,-0.5786415,0.4809155,0.2166510,-0.3291053,-0.5252405,-0.4780703,-1.3188585,0.4904822,-0.5234721,-1.6400329,-2.0174474,0.5317536,-0.0138196,-1.4635940,-1.8695927,-1.7503391;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x =
					ElementCreator::create(
						SpaceX->getDualSpace(),
						Xtemp.col(i));
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
	{
		const double power = Helpers::ConjugateValue(2.);
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						Xtemp.innerSize(), "lp", args);
		std::vector<SpaceElement_ptr_t> U =
				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
		const Mapping &J_p = *SpaceX->getDualSpace()->getDualityMapping();
		const BregmanProjectionFunctional d_p(
				*SpaceX->getDualSpace()->getNorm(),
				dynamic_cast<const PowerTypeDualityMapping &>(J_p),
				J_p.getPower(),U,alpha);
		Eigen::VectorXd expected(10);
		expected << 0.52896,0.41386,0.38071,1.51650,3.19301,0.56201,-0.10539,1.14000,1.10983,0.59547;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << -0.0364139,0.4636278,0.0099964,-0.6344602,-1.2532194,0.6526892,0.1775279,-0.3778562,-0.7041339,-0.3567879,-0.7014935,0.4905796,-0.5404139,-2.0966682,-3.6643320,0.9372717,-0.1213626,-1.5904568,-2.3005447,-1.4595545;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x =
					ElementCreator::create(
						SpaceX->getDualSpace(),
						Xtemp.col(i));
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
	{
		const double power = Helpers::ConjugateValue(4.);
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						Xtemp.innerSize(), "lp", args);
		std::vector<SpaceElement_ptr_t> U =
				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
		const Mapping &J_p = *SpaceX->getDualSpace()->getDualityMapping();
		const BregmanProjectionFunctional d_p(
				*SpaceX->getDualSpace()->getNorm(),
				dynamic_cast<const PowerTypeDualityMapping &>(J_p),
				J_p.getPower(),U,alpha);
		Eigen::VectorXd expected(10);
		expected << 0.51837,0.16386,0.23653,1.62694,7.79111,0.56620,-0.31918,0.91932,1.05715,0.39978;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << 0.013798,0.463720,0.013185,-1.414160,-6.837598,1.295795,0.115627,-0.512366,-1.342928,-0.183546,-0.571259,0.490796,-0.558629,-3.931688,-17.297786,2.455497,-0.291517,-1.940487,-3.839391,-1.044191;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x =
					ElementCreator::create(
						SpaceX->getDualSpace(),
						Xtemp.col(i));
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

void BregmanProjectionFunctionalUnitTest::inftyNorm()
{
	// due to large/small values we look at relative precision
	const double p = Helpers::ConjugateValue(std::numeric_limits<double>::infinity());
	Eigen::MatrixXd Utemp(2,2);
	Utemp <<   -0.37133, -0.95779, 0.40013, 0.93263;
	std::vector<double> alpha;
	alpha += 0.014897, -0.568410;
	Eigen::MatrixXd ttemp(2,10);
	ttemp << -0.104574,0.559468,0.766960,0.148054,-0.923356,0.041415,0.015722,-0.033294,-0.570823,-0.051950,-0.914503,0.166395,-0.342347,-0.728231,-0.888305,0.778963,0.730667,-0.825560,-0.306347,-0.579328;
	std::vector< std::vector<double> > t =
			getVectorOfVectors(ttemp);
	Eigen::MatrixXd Xtemp(2,10);
	Xtemp << 0.768586,-0.225475,0.470209,0.483057,-0.324979,-0.554573,0.030939,-0.352048,-0.103521,0.305334,-0.902129,-0.610977,0.396205,0.854957,0.562012,-0.658914,0.964865,-0.873672,0.717819,0.118345;
	{
		const double power = Helpers::ConjugateValue(1.1);
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						Xtemp.innerSize(), "lp", args);
		std::vector<SpaceElement_ptr_t> U =
				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
		const Mapping &J_infty = *SpaceX->getDualSpace()->getDualityMapping();
		const BregmanProjectionFunctional d_p(
				*SpaceX->getDualSpace()->getNorm(),
				dynamic_cast<const PowerTypeDualityMapping &>(J_infty),
				J_infty.getPower(),U,alpha);
		Eigen::VectorXd expected(10);
		expected << 0.62787,0.81287,0.56263,1.81008,2.18416,0.87617,0.23439,1.53412,1.30917,0.92277;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << -0.29146,0.41462,0.35594,-0.40109,-0.40850,0.42878,0.37504,-0.36182,-0.39367,-0.37006,-1.35863,0.36329,0.31127,-1.53800,-1.55528,0.39627,0.36054,-1.54011,-1.52070,-1.46569;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x =
					ElementCreator::create(
						SpaceX->getDualSpace(),
						Xtemp.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t[i],x) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t[i],x) )/expected(i) ) < tolerance);
			const std::vector<double> gval = d_p.gradient(t[i],x);
//			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << gval << ".\n";
			Eigen::VectorXd difference_vector = expected_gradient.col(i);
			for (int j=0;j<difference_vector.size();++j)
				difference_vector[j] -= gval[j];
			CPPUNIT_ASSERT( (difference_vector.norm())/(expected_gradient.col(i)).norm() < tolerance);
		}
	}
	{
		const double power = Helpers::ConjugateValue(2.);
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						Xtemp.innerSize(), "lp", args);
		std::vector<SpaceElement_ptr_t> U =
				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
		const Mapping &J_infty = *SpaceX->getDualSpace()->getDualityMapping();
		const BregmanProjectionFunctional d_p(
				*SpaceX->getDualSpace()->getNorm(),
				dynamic_cast<const PowerTypeDualityMapping &>(J_infty),
				J_infty.getPower(),U,alpha);
		Eigen::VectorXd expected(10);
		expected << 0.52893,0.40382,0.29723,1.50379,2.03986,0.54061,-0.14379,1.13591,0.92445,0.55934;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << -0.039371,0.411033,0.173492,-0.575250,-0.689305,0.575867,0.288415,-0.414030,-0.478035,-0.256964,-0.708385,0.354918,-0.159334,-1.943942,-2.209787,0.739118,0.137097,-1.674771,-1.717351,-1.202072;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x =
					ElementCreator::create(
						SpaceX->getDualSpace(),
						Xtemp.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t[i],x) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t[i],x) )/expected(i) ) < tolerance);
			const std::vector<double> gval = d_p.gradient(t[i],x);
//			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << gval << ".\n";
			Eigen::VectorXd difference_vector = expected_gradient.col(i);
			for (int j=0;j<difference_vector.size();++j)
				difference_vector[j] -= gval[j];
			CPPUNIT_ASSERT( (difference_vector.norm())/(expected_gradient.col(i)).norm() < tolerance);
		}
	}
	{
		const double power = Helpers::ConjugateValue(10.);
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						Xtemp.innerSize(), "lp", args);
		std::vector<SpaceElement_ptr_t> U =
				getVectorOfSpaceElements(SpaceX->getDualSpace(), Utemp);
		const Mapping &J_infty = *SpaceX->getDualSpace()->getDualityMapping();
		const BregmanProjectionFunctional d_p(
				*SpaceX->getDualSpace()->getNorm(),
				dynamic_cast<const PowerTypeDualityMapping &>(J_infty),
				J_infty.getPower(),U,alpha);
		Eigen::VectorXd expected(10);
		expected << 5.1825e-01,4.2115e-03,2.0604e-01,5.2869e+00,2.9000e+01,2.4914e+00,-4.1038e-01,8.9170e-01,9.7075e-01,3.3062e-01;
		Eigen::MatrixXd expected_gradient(2,10);
		expected_gradient << 1.4897e-02,3.8049e-01,1.5072e-02,-1.3199e+01,-6.4800e+01,8.3874e+00,3.8602e-02,-1.3447e+00,-2.6001e+00,2.5507e-03,-5.6841e-01,2.8373e-01,-5.6796e-01,-3.1368e+01,-1.5164e+02,1.8947e+01,-5.0726e-01,-4.0753e+00,-6.6636e+00,-5.9719e-01;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x =
					ElementCreator::create(
						SpaceX->getDualSpace(),
						Xtemp.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(t[i],x) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(t[i],x) )/expected(i) ) < 2e-4);
			const std::vector<double> gval = d_p.gradient(t[i],x);
//			std::cout << "# " << i << ": Expecting " << expected_gradient.col(i).transpose() << " and got " << gval << ".\n";
			Eigen::VectorXd difference_vector = expected_gradient.col(i);
			for (int j=0;j<difference_vector.size();++j)
				difference_vector[j] -= gval[j];
			CPPUNIT_ASSERT( (difference_vector.norm())/(expected_gradient.col(i)).norm() < tolerance);
		}
	}
}
