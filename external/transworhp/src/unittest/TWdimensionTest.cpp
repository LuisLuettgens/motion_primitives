#include <boost/test/unit_test.hpp>

#include "TWproblem.h"

using namespace tw;

BOOST_AUTO_TEST_CASE(TWdimensionTest) {

	TWdimension TWdim;

	// default Werte
	BOOST_CHECK_EQUAL(TWdim.n_dis, 2);
	BOOST_CHECK_EQUAL(TWdim.n_ode, 0);
	BOOST_CHECK_EQUAL(TWdim.n_ctrl, 0);
	BOOST_CHECK_EQUAL(TWdim.n_param, 0);
	BOOST_CHECK_EQUAL(TWdim.n_rand, 0);
	BOOST_CHECK_EQUAL(TWdim.n_neben, 0);
	BOOST_CHECK_EQUAL(TWdim.n_integral, 0);
	BOOST_CHECK_EQUAL(TWdim.n_zen, 0);

	TWdim.ID = "Raketenwagen";
	TWdim.n_dis = 11;
	TWdim.n_ode = 2;
	TWdim.n_ctrl = 3;
	TWdim.n_param = 5;
	TWdim.n_rand = 7;
	TWdim.n_neben = 9;
	TWdim.n_integral = 13;
	TWdim.n_zen = 15;

	BOOST_CHECK_EQUAL(TWdim.ID, "Raketenwagen");
	BOOST_CHECK_EQUAL(TWdim.n_dis, 11);
	BOOST_CHECK_EQUAL(TWdim.n_ode, 2);
	BOOST_CHECK_EQUAL(TWdim.n_ctrl, 3);
	BOOST_CHECK_EQUAL(TWdim.n_param, 5);
	BOOST_CHECK_EQUAL(TWdim.n_rand, 7);
	BOOST_CHECK_EQUAL(TWdim.n_neben, 9);
	BOOST_CHECK_EQUAL(TWdim.n_integral, 13);
	BOOST_CHECK_EQUAL(TWdim.n_zen, 15);

	// copy ctor
	TWdimension TWdim2 = TWdim;
	BOOST_CHECK_EQUAL(TWdim.ID, TWdim2.ID);
	BOOST_CHECK_EQUAL(TWdim.n_dis, TWdim2.n_dis);
	BOOST_CHECK_EQUAL(TWdim.n_ode, TWdim2.n_ode);
	BOOST_CHECK_EQUAL(TWdim.n_ctrl, TWdim2.n_ctrl);
	BOOST_CHECK_EQUAL(TWdim.n_param, TWdim2.n_param);
	BOOST_CHECK_EQUAL(TWdim.n_rand, TWdim2.n_rand);
	BOOST_CHECK_EQUAL(TWdim.n_neben, TWdim2.n_neben);
	BOOST_CHECK_EQUAL(TWdim.n_integral, TWdim2.n_integral);
	BOOST_CHECK_EQUAL(TWdim.n_zen, TWdim2.n_zen);
	BOOST_CHECK(TWdim.multinode == TWdim2.multinode);
	BOOST_CHECK(TWdim.BOXneben == TWdim2.BOXneben);

	// setMultinodes
	TWdim.setMultinodes(3);
	std::vector<int> m = {0,5,10};
	BOOST_CHECK(TWdim.multinode == m);
}
