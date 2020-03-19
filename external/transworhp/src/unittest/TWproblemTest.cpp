#include "TWfolder.h"
#include "TWproblem.h"

#include <boost/test/unit_test.hpp>

using namespace tw;

/* pruefen, ob alle Problemgroessen richtig uebernommen werden */
BOOST_AUTO_TEST_CASE(TWproblemTest1) {

	TWdimension TWdim;
	TWdim.ID = "Raketenwagen";
	TWdim.n_dis = 11;
	TWdim.n_ode = 2;
	TWdim.n_ctrl = 3;
	TWdim.n_param = 5;
	TWdim.n_rand = 7;
	TWdim.n_neben = 9;
	TWdim.n_integral = 13;
	TWdim.n_zen = 15;

	EmptyTransWorhpProblem ph(TWdim);

	BOOST_CHECK_EQUAL(ph.id, TWdim.ID);
	BOOST_CHECK_EQUAL(ph.n_dis, TWdim.n_dis);
	BOOST_CHECK_EQUAL(ph.n_ode, TWdim.n_ode);
	BOOST_CHECK_EQUAL(ph.n_ctrl, TWdim.n_ctrl);
	BOOST_CHECK_EQUAL(ph.n_param, TWdim.n_param);
	BOOST_CHECK_EQUAL(ph.n_rand, TWdim.n_rand);
	BOOST_CHECK_EQUAL(ph.n_neben, TWdim.n_neben);
	BOOST_CHECK_EQUAL(ph.n_integral, TWdim.n_integral);
	BOOST_CHECK_EQUAL(ph.n_zen, TWdim.n_zen);
}

/* Auswertung der Zielfunktion: bei EmptyTransWorhpProblem obj() = 0.0 */
BOOST_AUTO_TEST_CASE(TWproblemTest2) {

	TWdimension TWdim;
	TWdim.ID = "Raketenwagen";
	TWdim.n_dis = 11;
	TWdim.n_ode = 2;
	TWdim.n_ctrl = 3;
	TWdim.n_param = 5;
	TWdim.n_rand = 7;
	TWdim.n_neben = 9;
	TWdim.n_integral = 13;
	TWdim.n_zen = 15;

	EmptyTransWorhpProblem ph(TWdim);

	const auto obj = ph.obj();

	BOOST_CHECK_EQUAL(obj, 0.0);
}

/* Index Funktionen auf den WORHP Vektor testen - Trapez */
BOOST_AUTO_TEST_CASE(TWproblemTest3) {

	TWdimension TWdim;
	TWdim.ID = "Raketenwagen";
	TWdim.n_dis = 11;
	TWdim.n_ode = 2;
	TWdim.n_ctrl = 3;
	TWdim.n_param = 5;
	TWdim.n_rand = 7;
	TWdim.n_neben = 9;
	TWdim.n_integral = 13;
	TWdim.n_zen = 15;

	TWparameter twparam("example/transworhp.xml");
	twparam.solver = TransWORHP_type::fullDiscretization;
	twparam.twdiscretization = TWdiscretization(TWdiscretizationType::Trapez,1,0);

	TWfolder folder(&twparam, 0);

	EmptyTransWorhpProblem ph(TWdim);
	ph.setSolver(&twparam);

	folder.Add(&ph);
	folder.Init();

	for (int dis = 0; dis < TWdim.n_dis; dis++) {
		for (int ode = 0; ode < TWdim.n_ode; ode++) {
			const auto index = ph.x_index(dis, ode);
			BOOST_CHECK_EQUAL(index, (TWdim.n_ode + TWdim.n_ctrl) * dis + ode);
		}
		for (int ctrl = 0; ctrl < TWdim.n_ctrl; ctrl++) {
			const auto index = ph.u_index(dis, ctrl);
			BOOST_CHECK_EQUAL(index, (TWdim.n_ode + TWdim.n_ctrl) * dis + TWdim.n_ode + ctrl);
		}
	}
	for (int param = 0; param < TWdim.n_param; param++) {
		const auto index = ph.p_index(param);
		BOOST_CHECK_EQUAL(index, (TWdim.n_ode + TWdim.n_ctrl) * TWdim.n_dis + param);
	}
}

/* Index Funktionen auf den WORHP Vektor testen - Hermite-Simpson */
BOOST_AUTO_TEST_CASE(TWproblemTest4) {

	TWdimension TWdim;
	TWdim.ID = "Raketenwagen";
	TWdim.n_dis = 11;
	TWdim.n_ode = 2;
	TWdim.n_ctrl = 3;
	TWdim.n_param = 5;
	TWdim.n_rand = 7;
	TWdim.n_neben = 9;
	TWdim.n_integral = 13;
	TWdim.n_zen = 15;

	TWparameter twparam("example/transworhp.xml");
	twparam.solver = TransWORHP_type::fullDiscretization;//=0, multipleShooting=2, pseudospectral=4, pseudospectral_gauss};;
	twparam.twdiscretization = TWdiscretization(TWdiscretizationType::HermiteSimpson,2,1);

	TWfolder folder(&twparam, 0);

	EmptyTransWorhpProblem ph(TWdim);
	ph.setSolver(&twparam);

	folder.Add(&ph);
	folder.Init();

	for (int dis = 0; dis < TWdim.n_dis; dis++) {
		for (int ode = 0; ode < TWdim.n_ode; ode++) {
			const auto index = ph.x_index(dis, ode);
			BOOST_CHECK_EQUAL(index, (TWdim.n_ode + TWdim.n_ctrl) * 2 * dis + ode);
		}
		for (int ctrl = 0; ctrl < TWdim.n_ctrl; ctrl++) {
			const auto index = ph.u_index(dis, ctrl);
			BOOST_CHECK_EQUAL(index, (TWdim.n_ode + TWdim.n_ctrl) * 2 * dis + TWdim.n_ode + ctrl);
		}
	}
	for (int param = 0; param < TWdim.n_param; param++) {
		const auto index = ph.p_index(param);
		BOOST_CHECK_EQUAL(index, (TWdim.n_ode + TWdim.n_ctrl) * 2 * (TWdim.n_dis-1) + (TWdim.n_ode + TWdim.n_ctrl) + param);
	}
}

/* Index Funktionen auf den WORHP Vektor testen - explTW */
BOOST_AUTO_TEST_CASE(TWproblemTest5) {

	TWdimension TWdim;
	TWdim.ID = "Raketenwagen";
	TWdim.n_dis = 11;
	TWdim.n_ode = 2;
	TWdim.n_ctrl = 3;
	TWdim.n_param = 5;
	TWdim.n_rand = 7;
	TWdim.n_neben = 9;
	TWdim.n_integral = 13;
	TWdim.n_zen = 15;

	TWdim.multinode = {0, 2, 10};

	TWparameter twparam("example/transworhp.xml");
	twparam.solver = TransWORHP_type::multipleShooting;//=2, pseudospectral=4, pseudospectral_gauss};;

	TWfolder folder(&twparam, 0);

	EmptyTransWorhpProblem ph(TWdim);
	ph.setSolver(&twparam);

	folder.Add(&ph);
	folder.Init();

	// Index: 0 - 1. Mehrzielknoten
	int dis_index = 0;
	for (int ode = 0; ode < TWdim.n_ode; ode++) {
		const auto index = ph.x_index(dis_index, ode);
		BOOST_CHECK_EQUAL(index, ode);
	}
	for (int ctrl = 0; ctrl < TWdim.n_ctrl; ctrl++) {
		const auto index = ph.u_index(dis_index, ctrl);
		BOOST_CHECK_EQUAL(index, TWdim.n_ode + ctrl);
	}

	// Index: 1 - KEIN Mehrzielknoten
	dis_index = 1;
	for (int ctrl = 0; ctrl < TWdim.n_ctrl; ctrl++) {
		const auto index = ph.u_index(dis_index, ctrl);
		BOOST_CHECK_EQUAL(index, TWdim.n_ode + TWdim.n_ctrl + ctrl);
	}

	// Index: 2 - 2. Mehrzielknoten
	dis_index = 2;
	for (int ode = 0; ode < TWdim.n_ode; ode++) {
		const auto index = ph.x_index(dis_index, ode);
		BOOST_CHECK_EQUAL(index, TWdim.n_ode + TWdim.n_ctrl*2 + ode);
	}
	for (int ctrl = 0; ctrl < TWdim.n_ctrl; ctrl++) {
		const auto index = ph.u_index(dis_index, ctrl);
		BOOST_CHECK_EQUAL(index, TWdim.n_ode + TWdim.n_ctrl*2 + TWdim.n_ode + ctrl);
	}

	// Index: 3 - KEIN Mehrzielknoten
	dis_index = 3;
	for (int ctrl = 0; ctrl < TWdim.n_ctrl; ctrl++) {
		const auto index = ph.u_index(dis_index, ctrl);
		BOOST_CHECK_EQUAL(index, TWdim.n_ode + TWdim.n_ctrl*2 + TWdim.n_ode + TWdim.n_ctrl + ctrl);
	}

	// Index: 4 - KEIN Mehrzielknoten
	dis_index = 4;
	for (int ctrl = 0; ctrl < TWdim.n_ctrl; ctrl++) {
		const auto index = ph.u_index(dis_index, ctrl);
		BOOST_CHECK_EQUAL(index, TWdim.n_ode + TWdim.n_ctrl*2 + TWdim.n_ode + TWdim.n_ctrl*2 + ctrl);
	}

	// Index: 5 - KEIN Mehrzielknoten
	dis_index = 5;
	for (int ctrl = 0; ctrl < TWdim.n_ctrl; ctrl++) {
		const auto index = ph.u_index(dis_index, ctrl);
		BOOST_CHECK_EQUAL(index, TWdim.n_ode + TWdim.n_ctrl*2 + TWdim.n_ode + TWdim.n_ctrl*3 + ctrl);
	}

	// Index: 6 - KEIN Mehrzielknoten
	dis_index = 6;
	for (int ctrl = 0; ctrl < TWdim.n_ctrl; ctrl++) {
		const auto index = ph.u_index(dis_index, ctrl);
		BOOST_CHECK_EQUAL(index, TWdim.n_ode + TWdim.n_ctrl*2 + TWdim.n_ode + TWdim.n_ctrl*4 + ctrl);
	}

	// Index: 7 - KEIN Mehrzielknoten
	dis_index = 7;
	for (int ctrl = 0; ctrl < TWdim.n_ctrl; ctrl++) {
		const auto index = ph.u_index(dis_index, ctrl);
		BOOST_CHECK_EQUAL(index, TWdim.n_ode + TWdim.n_ctrl*2 + TWdim.n_ode + TWdim.n_ctrl*5 + ctrl);
	}

	// Index: 8 - KEIN Mehrzielknoten
	dis_index = 8;
	for (int ctrl = 0; ctrl < TWdim.n_ctrl; ctrl++) {
		const auto index = ph.u_index(dis_index, ctrl);
		BOOST_CHECK_EQUAL(index, TWdim.n_ode + TWdim.n_ctrl*2 + TWdim.n_ode + TWdim.n_ctrl*6 + ctrl);
	}

	// Index: 9 - KEIN Mehrzielknoten
	dis_index = 9;
	for (int ctrl = 0; ctrl < TWdim.n_ctrl; ctrl++) {
		const auto index = ph.u_index(dis_index, ctrl);
		BOOST_CHECK_EQUAL(index, TWdim.n_ode + TWdim.n_ctrl*2 + TWdim.n_ode + TWdim.n_ctrl*7 + ctrl);
	}

	// Index: 10 - 3. Mehrzielknoten
	dis_index = 10;
	for (int ode = 0; ode < TWdim.n_ode; ode++) {
		const auto index = ph.x_index(dis_index, ode);
		BOOST_CHECK_EQUAL(index, TWdim.n_ode + TWdim.n_ctrl*2 + TWdim.n_ode + TWdim.n_ctrl*8 + ode);
	}
	for (int ctrl = 0; ctrl < TWdim.n_ctrl; ctrl++) {
		const auto index = ph.u_index(dis_index, ctrl);
		BOOST_CHECK_EQUAL(index, TWdim.n_ode + TWdim.n_ctrl*2 + TWdim.n_ode + TWdim.n_ctrl*8 + TWdim.n_ode + ctrl);
	}

	for (int param = 0; param < TWdim.n_param; param++) {
		const auto index = ph.p_index(param);
		BOOST_CHECK_EQUAL(index, TWdim.n_ode + TWdim.n_ctrl*2 + TWdim.n_ode + TWdim.n_ctrl*8 + TWdim.n_ode + TWdim.n_ctrl + param);
	}
}
