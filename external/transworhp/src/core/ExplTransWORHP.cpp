#include "ExplTransWORHP.h"

#include "conversion.h"
#include "xmlio.h"

#include "butcher.h"
#include "TWfolder.h"
#include "TWproblem.h"
#include "twstatus.h"

#include "../base/vectortools.h"
#include "../base/exception.h"

#include <iomanip>
#include <sstream>
#include <vector>
#include <algorithm>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef WIN32
#include <windows.h>
#endif

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::stringstream;

namespace tw {

ExplTransWorhp::ExplTransWorhp(TransWorhpProblem* ph, TWparameter *twparam, vector<int> multi, vector<int> box_neben)
	: TransWorhp(ph,twparam),
	multinodes(std::move(multi)),
	pointer_p(nullptr),
	boxNB(std::move(box_neben)), n_boxNeben(boxNB.size()) {

	if (multinodes.size() > static_cast<size_t>(n_dis)) {
		MyStatus("Error", "too many multinodes: " + std::to_string(multinodes.size()) + ", n_dis=" +
		                      std::to_string(n_dis),
		         Status::ERR);
		exit(1);
	} else if (multinodes.size() == 0) {
		MyStatus("Error", "there are no multinodes!", Status::ERR);
		exit(1);
	}

	n_multinodes = static_cast<int>(multinodes.size());
	transworhp_type = TransWORHP_type::multipleShooting;

#ifdef _OPENMP
	if (twparameter->parallel) {
		tmp_ode_parallel = vector<double>(n_ode*n_dis);

		tmp_ctrl_1_parallel = vector<double>(n_ctrl*n_dis);
		tmp_ctrl_2_parallel = vector<double>(n_ctrl*n_dis);

		butcher_parallel.reserve(multinodes.size()-1);
		for (size_t i = 0; i < multinodes.size()-1; i++) {
			butcher_parallel.emplace_back(n_ode, n_ctrl, twparameter->abserr, twparameter->relerr);
		}

		if (viewer) {
			n_threads = std::min(static_cast<int>(multinodes.size())-2,omp_get_max_threads()-1);
		} else {
			n_threads = std::min(static_cast<int>(multinodes.size())-1,omp_get_max_threads());
		}

		if (n_threads == 0) {
			n_threads = 1;
		}

		string s("max. number of Threads: " + std::to_string(n_threads) );
		MyStatus("Parallel", s, Status::NORMAL);
	}
#endif
}



ExplTransWorhp::~ExplTransWorhp() {

#ifdef _OPENMP
	butcher_parallel.clear();
#endif
}




void ExplTransWorhp::GetSpace(int delta1, int delta2) {

	Delta1 = delta1;
	Delta2 = delta2;

	// check BOX-constraints
	const double auxInfty = 1e200; // da Infty an dieser Stelle noch nicht bekannt ist

	std::fill(tmp_ode_1.begin(),tmp_ode_1.end(),-auxInfty);
	std::fill(tmp_ode_2.begin(),tmp_ode_2.end(),+auxInfty);


	phase->x_boundary(tmp_ode_1.data(),tmp_ode_2.data());

	vector<int> aux;
	aux.reserve(n_dis);

	// suche Beschraenkungen -> boxNB
	for (int i = 0; i < n_ode; i++) {
		if (tmp_ode_1[i] != -auxInfty || tmp_ode_2[i] != +auxInfty) {
			aux.push_back(i);
		}
	}

	if (boxNB.size() > 0) {
		if (*std::max_element(boxNB.begin(),boxNB.end()) >= n_ode) {
			stringstream ss;
			ss << "Try to constrain state " << *std::max_element(boxNB.begin(),boxNB.end());
			MyStatus("BOXneben", ss.str(), Status::ERR);
		}
	}

	if (boxNB.size() != aux.size()) {

		stringstream ss;
		Status errorCode = Status::WARN;
		ss << "Total number of states: " << aux.size() << ", constrained: " << boxNB.size();

		if (boxNB.size() < aux.size()) {
			ss << " (not all states are constrained)";
		} else {
			errorCode = Status::ERR;
			MyStatus("BOXneben", "Error while specifying number of constrained states:",
			         Status::ERR);
		}

		MyStatus("BOXneben",ss.str(),errorCode);
	}

	DS_obj.useStructure(phase->obj_structure(DS_obj));
	DS_ode.useStructure(phase->ode_structure(DS_ode));
	DS_rand.useStructure(phase->rand_structure(DS_rand));
	DS_neben.useStructure(phase->neben_structure(DS_neben));
	DS_integral.useStructure(phase->integral_structure(DS_integral));

	/*//if (n_integral)
	{

	if (DS_integral.use_struct) {

	for (int k=0; k<n_integral; k++) {
	for (int i=0; i<twdiscretization->stuetzstellen(n_dis); i++) {
	for (int j=0; j<n_ode + n_ctrl; j++) {
	if (DS_integral.check(k,j)) {
	DS_obj(0, i * (n_ode+n_ctrl) + j);
	}
	}
	}

	for (int i=0; i<n_param; i++) {
	if (DS_integral.check(k,n_ode + n_ctrl + i)) {
	DS_obj(0, n_dis * (n_ode+n_ctrl) + i);
	}
	}
	}
	}

	}*/



	DS_integral.finish();

	DS_obj.finish();

	DS_ode.finish();
	while (DS_ode.ode_power(n_ode,n_ctrl)) {}

	DS_neben.finish();
	DS_neben.neben_power(DS_ode,n_ode);

	printStructure();
}


string ExplTransWorhp::type_G(int row) const {

	stringstream a;
	//TODO
	if (row < n_ode*(n_dis*n_multinodes-1)) {
		a << " COUNT " << std::setw(2) << row%n_ode << " ";
	}
	else if (row < n_ode*(n_dis*n_multinodes-1)+n_rand) {
		a << " RAND  " << std::setw(2) << (row-n_ode*(twdiscretization->stuetzstellen(n_dis)-1))%n_rand << " ";
	}
	//else if (row <  n_ode*(n_dis*n_multinodes-1)+n_rand+n_neben*(n_dis)) {
	//	a << " NB   " << std::setw(2) << (row-n_ode*(twdiscretization->stuetzstellen(n_dis)-1)-n_rand-n_neben)%n_neben << " ";
	//}

	return a.str();
}



double ExplTransWorhp::x(int dis, int ode) const {

	return pointer_x[dis][ode];

	/*vector<int>::const_iterator it = find(multinode.begin(), multinode.end(), dis);
	if (it==multinode.end()) {
		//cout << "Error: x_index: " << dis << "invalid, only "<< multinode << " allowed."<< endl;
		return block_x[n_ode * twdiscretization->punkte() * dis + ode];
	}
	else {
		int ct = it-multinode.begin();
		return X[(n_ode ) * ct + ode];
	}
	*/

}

double &ExplTransWorhp::setx(int dis, int ode) {

	//cout << "setx()" << endl;
	//return pointer_x[dis][ode];
	return block_x[n_ode * twdiscretization->punkte() * dis + ode];
}

double ExplTransWorhp::u(int dis, int ctrl) const {

	return pointer_u[dis][ctrl];
	//return X[n_ode*n_multinodes + n_ctrl * twdiscretization->punkte() * dis + ctrl];
}

double ExplTransWorhp::p(int param) const {

	return pointer_p[param];
	//return X[n_ode*n_multinodes + n_ctrl * twdiscretization->stuetzstellen(n_dis) + param];
}

// hier keine Unterscheidung wie bei eingeschobenen Stuetzstellen (H-S) noetig
double ExplTransWorhp::x__(int dis, int ode) const {
	return x(dis,ode);
}

double ExplTransWorhp::u__(int dis, int ctrl) const {
	return u(dis,ctrl);
}


int ExplTransWorhp::x_index(int dis, int ode) const {
	// x_index nur auf Stuetzstellen anwendbar!?!

	if (pointer_x.empty()) {
		// nachrechnen...
		int index_start = 0;
		int m_dis = 0; // zaehlt disk punkte
		for (int m_node = 0; m_node < n_multinodes-1; m_node++) {
			int area = n_ode + n_ctrl * (multinodes[m_node+1]-multinodes[m_node]);
			if (m_dis == dis) break;

			index_start += area;
			m_dis += (multinodes[m_node+1]-multinodes[m_node]);
		}

		return index_start + ode;
	}

	int index = &pointer_x[dis][ode] - X;
	if (index>=0 && index<n_var) {
		return index;
	}
	else {
		cout << "Error: x_index: " << dis << " invalid, only "<< multinodes << " allowed."<< endl;
		return 0;
	}
}

int ExplTransWorhp::x_index__(int dis, int ode) const {
	return x_index(dis,ode);
}

int ExplTransWorhp::u_index(int dis, int ctrl) const {

	if (pointer_u.empty()) {
		int m = 0;

		if (dis == 0) {
			if(std::binary_search(multinodes.begin(), multinodes.end(), 0)) {
				m += n_ode;
			}
			return m + ctrl;
		}
		else if (dis == n_dis-1) {
			return n_multinodes*n_ode+(n_dis-1)*n_ctrl+ctrl;
		}

		for (int i=0; i < dis+1; i++) {
			//Pruefen, ob Gitterpunkt Multiknoten ist
			if(std::binary_search(multinodes.begin(), multinodes.end(), i)) {
				m += n_ode;
			}
		}
		m += dis*n_ctrl;

		return m + ctrl;
	}

	int index = &pointer_u[dis][ctrl] - X;
	if (index>=0 && index<n_var) {
		return index;
	}
	else {
		cout << "Error: u_index: " << dis << " invalid."<< endl;
		return 0;
	}

	//return n_ode*n_multinodes + n_ctrl * twdiscretization->punkte() * dis  + ctrl;
}

int ExplTransWorhp::u_index__(int dis, int ode) const {
	return u_index(dis,ode);
}

int ExplTransWorhp::p_index(int param) const {

	if (!pointer_p) {
		return n_ode*n_multinodes + n_ctrl * twdiscretization->stuetzstellen(n_dis) + param;
	}

	int index = &pointer_p[param] - X;
	if (index>=0 && index<n_var) {
		return index;
	}
	else {
		cout << "Error: p_index: " << param << " invalid." << index << endl;
		return 0;
	}
}

int ExplTransWorhp::x_indexode(int ode) const {
	return ode;
}

int ExplTransWorhp::u_indexode(int ctrl) const {
	return n_ode + ctrl;
}

int ExplTransWorhp::p_indexode(int param) const {
	return n_ode + n_ctrl + param;
}


int ExplTransWorhp::dis_index(int dis) {

	if (T_ohneMulti.empty()) {
		T_ohneMulti.reserve(n_dis-n_multinodes);

		std::vector<int> auxT;
		auxT.reserve(n_dis);
		for (int i = 0; i < n_dis; i++) {
			auxT.push_back(i);
		}

		// Manuelles set_difference:
		/*
		auto it = multinode.begin();
		for (;it!=multinode.end();it++) {
			auto fit = find(auxT.begin(), auxT.end(), *it);
			if (fit!=auxT.end())
				auxT.erase(fit);
		}
		T_ohneMulti = auxT;
		*/
		//range-for
		for (int multi : multinodes) {
			auto fit = std::find(auxT.begin(), auxT.end(), multi);
			if (fit!=auxT.end())
				auxT.erase(fit);
		}
		T_ohneMulti = auxT;

		// Matthias K: set_difference macht bei mir unter VS2012 manchmal komische Sachen...
		//std::set_difference(auxT.begin(), auxT.end(), multinode.begin(), multinode.end(), T_ohneMulti.begin());
	}

	if (dis < n_dis-n_multinodes) {
		return T_ohneMulti[dis];
	} else {
		cout << "Error: dis_index: " << dis << " invalid." << endl;
  		return -1;
	}
  /*
	int aux = 0;

	for (int i = 0; i < n_dis; i++, aux++) {
		if(std::find(multinode.begin(), multinode.end(), i) != multinode.end()) {
			i++;
		}
		if (aux == dis) {
			return i;
		}
	}
	return -1;*/
}

bool ExplTransWorhp::isMultinode(int i) {

	if (i >= n_dis) {
		cout << "Error: isMultinode: " << i << " out of range." << endl;
		return -1;
	}

	//Werte fuellen, falls noch nicht geschehen
	if (isMulti.empty()) {
		isMulti.reserve(n_dis);

		for (int i = 0; i < n_dis; i++) {

			// binary_search, da die Multiknoten geordnet sind (bzw sein sollten!)
			if(std::binary_search(multinodes.begin(), multinodes.end(), i)) {
				isMulti.push_back(true);
			} else {
				isMulti.push_back(false);
			}
		}
	}

	return isMulti[i];
}


bool ExplTransWorhp::ode_diff2(DiffStructure&, double, const double*, const double*, const double*) {
	return false;
}


void ExplTransWorhp::Connect(const OptVar &o, const Params &p) {

	//cout << "TW::CONN" << endl;

	Infty = p.Infty;

	X = &o.X[Delta1];
	X_low = &o.XL[Delta1];
	X_upp = &o.XU[Delta1];


	G = &o.G[Delta2 + 0];
	G_low = &o.GL[Delta2 + 0];
	G_upp = &o.GU[Delta2 + 0];


	Lambda = &o.Lambda[Delta1];
	Mu = &o.Mu[Delta2];

	for (int i=0; i<n_var; i++) {
		Lambda[i] = 0.0;
	}
	for (int i=0; i<n_con; i++) {
		Mu[i] = 0.0;
	}

	std::fill(tmp_ode_1.begin(),tmp_ode_1.end(),0.0);
	std::fill(tmp_ctrl_1.begin(),tmp_ctrl_1.end(),0.0);

	// Speziell fuer Explizites Verfahren
	block_x = vector<double>(n_ode*twdiscretization->stuetzstellen(n_dis));
	block_u = vector<double>(n_ctrl*twdiscretization->stuetzstellen(n_dis));

	pointer_x = vector<double*>(twdiscretization->stuetzstellen(n_dis));
	pointer_u = vector<double*>(twdiscretization->stuetzstellen(n_dis));

	double *X0 = X;
	int dis = 0;
	for (int i = 0; i < twdiscretization->stuetzstellen(n_dis); i++) {
		if (dis < n_multinodes && multinodes[dis] == i) {
			pointer_x[i] = X0;
			X0 += n_ode;
			pointer_u[i] = X0;
			X0 += n_ctrl;
			dis++;
		}
		else {
			pointer_x[i] = &block_x[n_ode*i];
			pointer_u[i] = X0;
			X0 += n_ctrl;
		}
	}

	pointer_p = X0;


	ZEN = &o.P[0];
	phase->zen_init(ZEN);


	for (int dis = 0; dis < n_multinodes; dis++) {

		phase->x_init(tmp_ode_1.data(), multinodes[dis], twdiscretization->stuetzstellen(n_dis));

		for (int ode = 0; ode < n_ode; ode++) {
			pointer_x[multinodes[dis]][ode] = tmp_ode_1[ode];
		}
	}

	// graphic - initial state
	for (int dis = 0; dis < twdiscretization->stuetzstellen(n_dis); dis++) {
		phase->x_init(tmp_ode_1.data(), dis, twdiscretization->stuetzstellen(n_dis));
		for (int ode = 0; ode < n_ode; ode++) {
			block_x[dis * n_ode + ode] = tmp_ode_1[ode];
		}
	}

	for (int dis = 0; dis < twdiscretization->stuetzstellen(n_dis); dis++) {

		phase->u_init(tmp_ctrl_1.data(), dis, twdiscretization->stuetzstellen(n_dis));

		for (int ctrl = 0; ctrl < n_ctrl; ctrl++) {
			pointer_u[dis][ctrl] = tmp_ctrl_1[ctrl];
			// graphic - initial control
			block_u[dis * n_ctrl + ctrl] = tmp_ctrl_1[ctrl];
		}
	}

	phase->p_init(pointer_p);

	phase->init();
}


void ExplTransWorhp::Boundary() {

	std::fill(tmp_ode_1.begin(),tmp_ode_1.end(),-Infty);
	std::fill(tmp_ode_2.begin(),tmp_ode_2.end(),+Infty);

	std::fill(tmp_ctrl_1.begin(),tmp_ctrl_1.end(),-Infty);
	std::fill(tmp_ctrl_2.begin(),tmp_ctrl_2.end(),+Infty);


	phase->x_boundary(tmp_ode_1.data(),tmp_ode_2.data());

	phase->u_boundary(tmp_ctrl_1.data(),tmp_ctrl_2.data());

	for (int dis=0; dis<n_multinodes; dis++) {
		int index = pointer_x[multinodes[dis]] - X;
		for (int ode=0; ode<n_ode; ode++) {
			X_low[index + ode] = tmp_ode_1[ode];
			X_upp[index + ode] = tmp_ode_2[ode];
		}
	}

	for (int dis=0; dis<twdiscretization->stuetzstellen(n_dis); dis++) {
		int index = pointer_u[dis] - X;
		for (int ctrl=0; ctrl<n_ctrl; ctrl++) {
			X_low[index + ctrl] = tmp_ctrl_1[ctrl];
			X_upp[index + ctrl] = tmp_ctrl_2[ctrl];
		}
	}

	int index = pointer_p - X;
	for (int i=0; i<n_param; i++) {
		X_low[index+i] = -Infty;
		X_upp[index+i] = +Infty;
	}
	phase->p_boundary( &X_low[index], &X_upp[index] );

	phase->var_boundary(X_low, X_upp);

	// damit die Box-Schranken durch den Startwert nicht verletzt sind
	boxConToInitGuess();

	for (int i=0; i<(n_multinodes-1) * n_ode; i++) {
		G_low[i] = 0.0;
		G_upp[i] = 0.0;
	}

	std::fill(tmp_rand_1.begin(),tmp_rand_1.end(),0.0);
	std::fill(tmp_rand_2.begin(),tmp_rand_2.end(),0.0);

	if (n_rand) phase->rand_boundary(tmp_rand_1.data(),tmp_rand_2.data());

	for (int i=0; i<n_rand; i++) {
		G_low[ (n_multinodes-1) * n_ode + i] = tmp_rand_1[i];
		G_upp[ (n_multinodes-1) * n_ode + i] = tmp_rand_2[i];
	}

	std::fill(tmp_neben_1.begin(),tmp_neben_1.end(),-Infty);
	std::fill(tmp_neben_2.begin(),tmp_neben_2.end(),+Infty);

	if (n_neben) phase->neben_boundary(tmp_neben_1.data(),tmp_neben_2.data());

	for (int i=0; i<twdiscretization->stuetzstellen(n_dis); i++) {
		for (int j=0; j<n_neben; j++) {
			G_low[ (n_multinodes-1) * n_ode + n_rand + i*n_neben + j] = tmp_neben_1[j];
			G_upp[ (n_multinodes-1) * n_ode + n_rand + i*n_neben + j] = tmp_neben_2[j];
		}
	}

	for (int i=0; i<(twdiscretization->stuetzstellen(n_dis)-n_multinodes); i++) {
		int aux = (n_multinodes-1) * n_ode + n_rand + twdiscretization->stuetzstellen(n_dis)*n_neben;
		aux += i*n_boxNeben;
		for (int j=0; j<n_boxNeben; j++) {
			G_low[aux + j] = tmp_ode_1[boxNB[j]];
			G_upp[aux + j] = tmp_ode_2[boxNB[j]];
		}
	}

}


void ExplTransWorhp::Constraints2(double *GG, int DGflag) {

	if (DGflag == 0) {
		#ifndef _OPENMP
		IntegrateStates();
		#else
		if (twparameter->parallel) {
			IntegrateStates_parallel();
		} else {
			IntegrateStates();
		}
		#endif
	}

	for (int dis=0; dis<n_multinodes-1; dis++) {

		double *G1 = &GG[dis*n_ode];
		Continuous(G1, dis+1);
	}


	if (n_rand) {
		twcount_calls.rand.f++;
		phase->rand( &GG[n_ode*(n_multinodes-1)] );
	}

	if (n_neben) {
		for (int dis=0; dis<n_dis; dis++) {
			double *G1 = &GG[n_ode*(n_multinodes-1) + n_rand + dis*n_neben];

			double *X1 = pointer_x[dis];
			double *U1 = pointer_u[dis];
			double *P =  pointer_p;
			twcount_calls.neben.f++;

			phase->neben(G1,T[dis],X1,U1,P);
		}
	}

	if (n_boxNeben) {
		for (int dis=0; dis<(n_dis-n_multinodes); dis++) {
			double *G1 = &GG[n_ode*(n_multinodes-1) + n_rand + n_dis*n_neben + dis*n_boxNeben];

			double *X1 = pointer_x[dis_index(dis)];

			boxNeben(G1,X1);
		}
	}

}



int ExplTransWorhp::DG_structure(const TWfolder *f, WorhpMatrix *DG, int offset) {

	DG_start = offset;
	int ind = offset;

	for (int k=0; k<n_var; k++) {

		int dis = 0;
		int index_start = 0;

		for (;dis<n_multinodes-1;dis++) {
			int area = n_ode + n_ctrl * (multinodes[dis+1]-multinodes[dis]);
			if (k < index_start + area) break;
			index_start += area;
		}
		// dis = index des vorherigen Multiknotens zu k
		// index_start = index des WORHP-Variablen-Beginns des Multiknotens zu k

		int k_rel = k - index_start;
		// k_rel: worhp_o.X[k = index_start + k_rel]=x(multinode[dis],k_rel), wenn k_rel\in[0..n_ode-1]

		int u_ind=0;
		if (n_ctrl) u_ind = (k_rel - n_ode) / n_ctrl;
		// u_ind: worhp_o.X[k] = u(multinode[dis] + u_ind,(k_rel-n_ode) % n_ctrl)    wenn >=0


		//cerr << "KK " << k << " " << dis << " " << index_start << " mnode" << multinode[dis] << " " << k_rel << " u" << u_ind << endl;


		//cerr << "n_multi:" << n_multinodes << endl;

		for (int l=0; l<n_con; l++) {

			int add = 0;

			if (l<(n_multinodes-1)*n_ode) { // cont-block

				int c1 = l/n_ode;
				int c2 = l%n_ode;

				if (k_rel>=0 && k_rel < n_ode) { // x

					int x1 = dis;
					int x2 = k_rel;

					if (c1==x1) {
						int block_index=0, block_eq=0;
						double val, mult;

						if (DS_ode.odecheck(c2,x2, block_index, block_eq, val, mult, twdiscretization->type))
							add = 1;

					}
					if (c1+1==x1) { // Zielknoten
						if (c2==x2) add = 1;
					}

				}
				else if (k_rel>=0 && k < n_multinodes*n_ode + n_dis*n_ctrl) { // u

					int x1 = dis;
					int x2 = (k_rel-n_ode)%n_ctrl + n_ode;


					if (c1==x1) {
						int block_index=0, block_eq=0;
						double val, mult;
						if (DS_ode.odecheck(c2,x2, block_index, block_eq, val, mult, twdiscretization->type))
							add = 1;
					}

					if (c1+1==x1 && twparameter->linInter) {
						if (u_ind==0) {
							int block_index=0, block_eq=0;
							double val, mult;
							if (DS_ode.odecheck(c2,x2, block_index, block_eq, val, mult, twdiscretization->type))
								add = 1;
						}
					}
					//add=1;
				}
				else { // p

					int x2 = (k-n_multinodes*n_ode - n_dis*n_ctrl) + n_ctrl + n_ode;
					int block_index=0, block_eq=0;
					double val, mult;
					if (DS_ode.odecheck(c2,x2, block_index, block_eq, val, mult, twdiscretization->type))
						add = 1;

							//add = 1;
				}


			}
			else if (l<(n_multinodes-1)*n_ode + n_rand) { // Rand-Bedingungen

				if (DS_rand.check(l-(n_multinodes-1)*n_ode,k))
					add=1;
			}

			else if (l<(n_multinodes-1)*n_ode + n_rand + twdiscretization->stuetzstellen(n_dis) * n_neben) { // Neben-Bedingungen

				int c1 = (l - (n_multinodes-1)*n_ode - n_rand) / n_neben;
				int c2 = (l - (n_multinodes-1)*n_ode - n_rand) % n_neben;

				int x1 = dis;
				int x2 = 0;
				if (n_ctrl) x2 = (k_rel-n_ode)%n_ctrl + n_ode;

				int u_multinode = multinodes[x1];
				int u_multinode2 = -1;

				if (x1 < static_cast<int>(multinodes.size())-1) {
					u_multinode2 = multinodes[x1+1];
				}

				if (c1==u_multinode) { // Mehrzielknoten
					if (k_rel < n_ode+n_ctrl) {
						if (DS_neben.check(c2,k_rel)) {
							add=1;
						}
					}
				}

				if (c1>u_multinode && (u_multinode==-1 || c1<u_multinode2)) {

					if (k_rel < n_ode + (c1-u_multinode+1) * n_ctrl) {
						//cout << c2 << " " << k_rel << " " << x2 << endl;
						if (k_rel < n_ode) {
							if (DS_neben.check(c2,k_rel)) {
								add=1;
							}
						} else {
							if (DS_neben.check(c2,x2)) {
								add=1;
							}
						}
					}
				}

				if (k >= n_multinodes*n_ode + n_dis*n_ctrl) {// freie Parameter
					add=1;
				}

			//add=1;
			}

			else {// BOXneben

				//c1 = dis. Punkt
				int c1 = (l - (n_multinodes-1)*n_ode - n_rand - twdiscretization->stuetzstellen(n_dis)*n_neben) / n_boxNeben;
				//c2 = Box_NB
				int c2 = (l - (n_multinodes-1)*n_ode - n_rand - twdiscretization->stuetzstellen(n_dis)*n_neben) % n_boxNeben;
				//cerr << c1 << " " << dis_index(c1) << " " << c2 << endl;

				// Trafo auf Gitterpunkte ohne Multiknoten
				c1 = dis_index(c1);
				// Trafo auf Beschraenkung
				c2 = boxNB[c2];

				int u_multinode = multinodes[dis];
				int u_multinode2 = -1;

				if (dis < static_cast<int>(multinodes.size())-1) {
					u_multinode2 = multinodes[dis+1];
				}

				if (k_rel>=0 && k_rel < n_ode) { // x

					//int x1 = dis_index(dis);
					int x2 = k_rel;

					if (c1>u_multinode && (u_multinode2==-1 || c1<u_multinode2)) {
						int block_index=0, block_eq=0;
						double val, mult;

						if (DS_ode.odecheck(c2,x2, block_index, block_eq, val, mult, twdiscretization->type)) {
							add = 1;
						}
					}

				}
				else if (k_rel>=0 && k_rel < n_ode+n_ctrl + (c1-u_multinode) * n_ctrl) {

					//int x1 = dis_index(dis);
					int x2 = (k_rel-n_ode)%n_ctrl + n_ode;

					if (c1>u_multinode && (u_multinode2==-1 || c1<u_multinode2)) {
						int block_index=0, block_eq=0;
						double val, mult;

						if (DS_ode.odecheck(c2,x2, block_index, block_eq, val, mult, twdiscretization->type)) {
							add = 1;
						}
					}

					//add=1;
				}

				if (k >= n_multinodes*n_ode + n_dis*n_ctrl) { // p
					add = 1;
				}
			//add= 1;
			}

			//if (f) f->DG_structure(ind,colindex + Delta1, DG);

			//add = 1;

			//if (DS_obj.check(l,k))
			if (add) {
				if (DG) {
					DG->col[ind] = k+1+Delta1;
					DG->row[ind] = l+1+Delta2;
					DG->val[ind] = 0;
				}
				ind++;
			}
		}
		if(f) { //Struktur aus TWfolder
			//cout << "Struktur aus TWfolder" << endl;
			f->DG_structure(ind, k + Delta1 + 1, DG);
		}
	}

	DG_nnz = ind - DG_start;

	//cout << "DG_STR" << DG_nnz << endl;
	return ind;
}


void ExplTransWorhp::DG_calculate(TWfolder *f, WorhpMatrix &DG) {

	int used = -1; //damit die NB nur einmal pro Obt.variable ausgewertet werden muessen

	for (int i = DG_start + DG_nnz-1; i >= DG_start; --i) {

		const int var = DG.col[i]-1-Delta1; //-1 wegen Fortran-Index
		const int con = DG.row[i]-1-Delta2;

		int start = typeOfX.at(var).dis;
		int end = n_dis-1;

		// wegen linearer Interpolation: vorheriges Intervall auch integrieren
		if (start > 0) {
			start -= 1;
		}

		// wenn es nicht das letzte Intervall ist, nur bis zum naechsten Multiknoten integrieren
		if (typeOfX.at(var).dis != n_dis-1) {
			if (n_multinodes > 1) {
				end = multinodes[typeOfX.at(var).multinodeL+1];
			}
		}

		// wenn ein Zustand gestoert wird, reicht es ab dort zu integrieren
		if (typeOfX.at(var).type == 0) {
			start = typeOfX.at(var).dis;
		}

		// fuer freie Parameter muss alles integriert werden
		else if (typeOfX.at(var).type == 2) {
			start = 0;
			end = n_dis-1;
		}


		if (con < n_con) {

			if (used != var) {

				// zentraler Differenzenquotient
				const double tmp = X[var];
				X[var] += eps;

				IntegrateStates2(start,end,1);
				Constraints2(tmp_gg_1.data(),1);

				X[var] = tmp - eps;
				IntegrateStates2(start,end,1);
				Constraints2(tmp_gg_2.data(),1);

				X[var] = tmp;
			}

			DG.val[i] = (tmp_gg_1[con]-tmp_gg_2[con])/(2*eps);
			used = var;
		} else { //Beschraenkungen aus dem TWfolder
			if (f) {
				int ind = i-DG_start;
				f->DG_calculate(DG_start,ind,DG.col[i],DG);
			}
		}
	}
}


int ExplTransWorhp::HM_structure_ohne_Diag(int hessianstructure , WorhpMatrix& /*DF*/, WorhpMatrix& /*DG*/, WorhpMatrix* HM, int offset) {
	//cout << "Expl::HM_structure_ohne_Diag" << HM_start << endl;

	HM_start = offset;

	int ind = offset;

	if (hessianstructure==0) { // Diagonal

	}

	else if (hessianstructure==1) { // Full

		for (int y = 0; y < n_var; ++y) {
			for (int x = y+1; x < n_var; ++x) {
				if (HM) {
					HM->row[ind] = x + 1 + Delta1;
					HM->col[ind] = y + 1 + Delta1;
				}
				ind++;
			}
		}
	}

	else if (hessianstructure==2) { // Blockode2, erweitert auf benachbarte Bloecke

		cout << "hessianstructure 2 " << endl;

		/*
		for (int y = 0; y < n_var; ++y) {
			for (int x = y+1; x < n_var; ++x) {

				int use = 0;

				int xblock = x / ((n_ode+n_ctrl));
				int yblock = y / ((n_ode+n_ctrl));

				if (abs(xblock-yblock)<2) use=1; // Benachbarte Bloecke verbinden

				if (x >= (n_ode+n_ctrl)*twdiscretization->stuetzstellen(n_dis) ) use=1; // mit freien Par. verbinden

				if (use) {
					if (HM) HM->row[ind] = x + 1 + Delta1;
					if (HM) HM->col[ind] = y + 1 + Delta1;
					ind++;
				}
			}
		}
		*/
	}

	else if (hessianstructure==3) { // Blockode

		//Multiknoten
		int disY = 0;
		//Blockgroesse
		int blocksizeY = 0;
		//Anfang des Blocks (in Y-Richtung)
		int index_startY = 0;

		//Groesse der Bloecke
		/*
		int blocksize[n_multinodes-1];
		for (int i = 0; i < n_multinodes-1; i++) {
			blocksize[i] = 2*n_ode + n_ctrl*(multinode[i+1]-multinode[i]+1);
			cout << "Block" << i << " " << blocksize[i] << endl;
		}
		*/

		for (int y = 0; y < n_var-1; y++) {

			for (;disY<n_multinodes-1;disY++) {
				//Blockgroesse minus Ueberschneidung
				blocksizeY = n_ode + n_ctrl * (multinodes[disY+1]-multinodes[disY]);
				if (y < index_startY + blocksizeY) break;
				index_startY += blocksizeY;
			}

			if (disY != n_multinodes-1) { //alle Bloecke vor dem letzten Multiknoten

				for (int x = y+1; x < index_startY+blocksizeY; x++) {

					if (x > n_var-n_param-1) break;

					if (HM) {
						HM->row[ind] = x + 1 + Delta1;
						HM->col[ind] = y + 1 + Delta1;
					}
					ind++;
				}

				// zwischen den for werden Zustaende weggelassen

				for (int x = index_startY+blocksizeY+n_ode; x < index_startY+(blocksizeY+n_ode+n_ctrl); x++) {

					if (x > n_var-n_param-1) break;

					if (HM) {
						HM->row[ind] = x + 1 + Delta1;
						HM->col[ind] = y + 1 + Delta1;
					}
					ind++;
				}
			} else { //letzter Multiknoten
				for (int x = y+1; x < n_var-n_param; x++) {

					if (HM) {
						HM->row[ind] = x + 1 + Delta1;
						HM->col[ind] = y + 1 + Delta1;
					}
					ind++;
				}
			}
			//freie Parameter
			for (int x = n_var-n_param; x < n_var; x++) {
				if (x > y) {
					if (HM) {
						HM->row[ind] = x + 1 + Delta1;
						HM->col[ind] = y + 1 + Delta1;
					}
					ind++;
				}
			}
		}
	}

	HM_nnz_ohne_diag = ind - HM_start;

	return ind;
}


void ExplTransWorhp::HM_calculate1(WorhpMatrix& /*DF*/, WorhpMatrix& /*DG*/, WorhpMatrix& /*HM*/, double /*ScaleObj*/, double* /*Mu*/) {
	//hier nicht moeglich
}


void ExplTransWorhp::Lagrange() {

	// fuer Trapezregel:
	for (int k=0; k<n_integral; k++) {
		lagrange_integral[k] = 0;
	}
	return;

	/* Bolza-Form für explTW implementieren
	double *P = &X[(n_ode+n_ctrl)* twdiscretization->stuetzstellen(n_dis)];

	if (twdiscretization->type==0) {
		for (int dis=0; dis<n_dis-1; dis++) {

			double *X1 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis];
			double *U1 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis + n_ode];

			double *X2 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1)];
			double *U2 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1) + n_ode];

			integral(tmp_integral_1.data(), T[dis], X1, U1, P);
			integral(tmp_integral_2.data(), T[dis+1], X2, U2, P);

			for (int k=0; k<n_integral; k++) {
				lagrange_integral[n_integral*(dis+1) + k ] = lagrange_integral[n_integral*(dis) + k ]
				+ (tmp_integral_1[k]+tmp_integral_2[k])/2 * (T[dis+1]-T[dis]);
			}
		}

	}

	if (twdiscretization->type==1) {
		for (int dis=0; dis<n_dis-1; dis++) { // nur Intervalle abgehen!

			double *X1 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis];
			double *U1 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*dis + n_ode];

			double *X2 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1)];
			double *U2 = &X[(n_ode+n_ctrl)* twdiscretization->punkte()*(dis+1) + n_ode];

			double *X12 = &X[(n_ode+n_ctrl)* (twdiscretization->punkte()*dis+1)];
			double *U12 = &X[(n_ode+n_ctrl)* (twdiscretization->punkte()*dis+1) + n_ode];

			double t1 = T[dis];
			double t2 = T[dis+1];
			double t12 = (t1+t2)/2;

			integral(tmp_integral_1.data(), t1, X1, U1, P);
			integral(tmp_integral_12, t12, X12, U12, P);
			integral(tmp_integral_2, t2, X2, U2, P);

			for (int k=0; k<n_integral; k++) {
				lagrange_integral[n_integral*(2*dis+2) + k] = lagrange_integral[n_integral*(2*dis)+k]
				+ (tmp_integral_1[k] + 4*tmp_integral_12[k] + tmp_integral_2[k])/6 * (t2-t1);

				lagrange_integral[n_integral*(2*dis+1) + k] = lagrange_integral[n_integral*(2*dis)+k]
				+ (5*tmp_integral_1[k] + 8*tmp_integral_12[k] - tmp_integral_2[k])/24 * (t2-t1);

			}
		}

	}
	*/


}


int ExplTransWorhp::Integrate(int /*btableau*/) {
	twdiscretization->type = TWdiscretizationType::MultipleShooting;
	// Integrate!
	for (int i=0;i<n_ode;i++) {
		block_x[i] = X[i];
	}

	return IntegrateInitial();
	//return TransWorhp::Integrate(btableau); // Methode aus TransWORHP (normal)

}

void ExplTransWorhp::Continuous(double* ddx, int dis) {

	//cout << "CONT " << pointer_x[multinode[dis]]-X <<endl;

	for (int i=0; i<n_ode; i++) {
		// pointer_x[dis] zeigt hoff. auf Opt-Variablen!
		ddx[i] = block_x[multinodes[dis] * n_ode + i] - pointer_x[multinodes[dis]][i];
	}
}

void ExplTransWorhp::IntegrateStates() {

	// fuer Grafik...
	for (int i = 0; i < n_ode; i++) {
		block_x[i] = X[i];
	}

	if (!schrittweiten.empty()) {
		const size_t auxSize = schrittweiten.size(); // alte Groesse als Richtwert fuer neue
		schrittweiten.clear();
		schrittweiten.reserve(auxSize);
	}

	// Anzahl der Integrationsschritte
	int steps = 0;
	for (int i = 0; i < n_dis-1; i++) {
		steps += integrate_schrittweite(i);
	}

	twcount_calls.integrate.f += steps;

	// fuer Grafik...
	for (int i = 0; i < n_dis; i++) {
		for (int j=0;j<n_ctrl;j++) {
			block_u[i*n_ctrl + j] = pointer_u[i][j];
		}
	}

}

void ExplTransWorhp::IntegrateStates2(int start, int end, int DGflag) {

	if (start < 0 || end > n_dis-1) {
		cout << "Error: IntegrateStates2: invalid input ["<< start << "," << end << "]." << endl;
	}

	// Anzahl der Integrationsschritte
	int steps = 0;
	for (int i=start; i<end; i++) {
		steps += integrate_schrittweite(i,DGflag);
	}

	twcount_calls.integrate.f += steps;
}

#ifdef _OPENMP
void ExplTransWorhp::IntegrateStates_parallel() {

	// fuer Grafik...
	for (int i = 0; i < n_ode; i++) {
		block_x[i] = X[i];
	}

	if (schrittweiten.empty()) {
		schrittweiten.resize(n_dis);
	}

	// Anzahl der Integrationsschritte
	int steps = 0;
	#pragma omp parallel for schedule(static) reduction(+:steps) num_threads(n_threads)
	for (int m = 1; m < n_multinodes; ++m) {
		for (int i = multinodes[m-1]; i < multinodes[m]; ++i) {
			steps += integrate_schrittweite_parallel(m-1,i);
		}
	}

	twcount_calls.integrate.f += steps;

	// fuer Grafik...
	for (int i = 0; i < n_dis; i++) {
		for (int j=0;j<n_ctrl;j++) {
			block_u[i*n_ctrl + j] = pointer_u[i][j];
		}
	}
}
#endif

int ExplTransWorhp::IntegrateInitial() {

	// fuer Grafik...
	for (int i=0;i<n_ode;i++) {
		block_x[i] = X[i];
	}

	int steps = 0; // Integrationsschritte

	for (int i=0; i<n_dis-1; i++) {

		steps += integrate_schrittweite(i);

		// setzen der Zustaende an den Multiknoten
		for (int ode=0;ode<n_ode;ode++) {
			pointer_x[i+1][ode] = setx(i+1,ode);
		}
	}

	// fuer Grafik...
	for (int i=0; i<n_dis; i++) {
		for (int j=0;j<n_ctrl;j++) {
			block_u[i*n_ctrl + j] = pointer_u[i][j];
		}
	}

	return steps;
}



double ExplTransWorhp::Objective(double ScaleObj) {

	IntegrateStates();

	twcount_calls.obj.f++;

	Lagrange();
	double a = phase->obj();
	for (int i=0; i<n_integral; i++)
		a += lagrange_weight[i] * lagrange_integral[n_integral*(twdiscretization->stuetzstellen(n_dis)-1)+i];
	//cout << a << endl;
	return a * ScaleObj;
}



int ExplTransWorhp::DF_structure(WorhpMatrix *DF, int offset) {

	DF_start = offset;
	int ind = offset;
	for (int k=0; k<n_var; k++) {

		if (DS_obj.check(0,k))
		{
			if (DF) {
				DF->row[ind] = k+1 + Delta1;
				DF->val[ind] = 0;
			}
			ind++;
		}
	}

	DF_nnz = ind - DF_start;

	return ind;
}


void ExplTransWorhp::DF_calculate(WorhpMatrix &DF, double ScaleObj) {

	for (int k=0; k<DF.nnz; k++) {

		int index = DF.row[k]-1;
		double v = X[index];
		X[index] += eps;
		IntegrateStates();
		twcount_calls.obj.f+=2;
		Lagrange();
		double f1 = phase->obj();
		for (int i=0; i<n_integral; i++)
			f1 += lagrange_weight[i] * lagrange_integral[n_integral*(twdiscretization->stuetzstellen(n_dis)-1)+i];

		X[index] = v- eps;
		IntegrateStates();
		Lagrange();
		double f2 = phase->obj();
		for (int i=0; i<n_integral; i++)
			f2 += lagrange_weight[i] * lagrange_integral[n_integral*(twdiscretization->stuetzstellen(n_dis)-1)+i];


		DF.val[k] = ScaleObj * (f1-f2)/(2*eps);

		X[index] = v;

	}
}

void ExplTransWorhp::boxNeben(double *c, const double *x) {

	for (int i = 0; i < n_boxNeben; i++) {
		c[i] = x[boxNB[i]];
	}
}

void ExplTransWorhp::boxNeben_boundary(double *c_low, double *c_upp) {

	for (int i = 0; i < n_boxNeben; i++) {
		c_low[i] = X_low[boxNB[i]];
		c_upp[i] = X_upp[boxNB[i]];
	}
}



int ExplTransWorhp::DoubleFrom(TransWorhp */*ph*/) {

	//muss hier neu implementiert werden,
	//da x_index, etc nicht ueberall definiert ist, wg Multiknoten

	cout << "DoubleFrom" << endl;
	exit(1);
	return -1;
}


void ExplTransWorhp::fromMATLAB_impl(std::ifstream &stream) {

	string line;

	/// Header
	getline(stream, line);
	//cout << ":" << line << endl;

	/// Parameter p
	getline(stream, line);
	//cout << ":" << line << endl;
	string line2 = string(line,2);
	vector<double> v = ToDoubleArray(line2);
	for (int i = 0; i < n_param; i++) {
		X[p_index(i)] = v[i];
	}

	/// Data
	getline(stream, line);
	//cout << ":" << line << endl;

	v = ToDoubleArray(line);
	if (phase->freetime) v[0] = v[0]/p(0);

	vector<double> v_last = v;
	//cout << ":" << v_last << " " << n_dis << endl;
	for (int i = 0; i < n_dis; i++) {

		while (T[i] > v[0]) {
			//cout << i << "   " << T[i] << " " << v[0] << endl;
			getline(stream, line);

			if (!stream) {
				v_last = v;
				break;
			} else {
				v_last = v;

				v = ToDoubleArray(line);
				if (phase->freetime) v[0] = v[0]/p(0);
			}
		}

		double vv = 1;

		if (v[0] != v_last[0]) {
			vv = (T[i] - v_last[0]) / (v[0] - v_last[0]);
		}
		//cout << i << "   vv " << vv << endl;

		for (int j = 0; j < n_ode; j++) {

			if (isMultinode(i)) {
				X[x_index(i,j)] = v_last[j+1] + (v[j+1]-v_last[j+1]) * vv;
			}

			// fuer Grafik...
			block_x[i*n_ode + j] = v_last[j+1] + (v[j+1]-v_last[j+1]) * vv;
		}

		for (int j = 0; j < n_ctrl; j++) {
			X[u_index(i,j)] = v_last[j+n_ode+1] + (v[j+n_ode+1]-v_last[j+n_ode+1]) * vv;

			// fuer Grafik...
			block_u[i*n_ctrl + j] = X[u_index(i,j)];
		}
	}

#ifdef TRANSWORHP_GRAPHICS
	if (viewer) viewer->autoScale();
#endif
}


void ExplTransWorhp::GetState(double *xx, double t) {

	int i=0;
	for (; i<n_dis; i++) {
		if (t < T[i]) break;
	}

	for (int j=0; j<n_ode; j++) {
		if (i>0 && i<n_dis)
			xx[j] = x(i-1,j) +  (x(i,j) - x(i-1,j)) / (T[i]-T[i-1])*(t-T[i-1]);
		if (i>=n_dis)
			xx[j] = x(n_dis-1,j);
	}

	//cout << t << " " << 0 << " " << xx[0] << endl;
}


void ExplTransWorhp::GetControl(double *uu, double t) {

	int i=0;
	for (; i<n_dis; i++) {
		if (t < T[i]) break;
	}

	for (int j=0; j<n_ctrl; j++) {
		if (i>0 && i<n_dis)
			uu[j] = u(i-1,j) +  (u(i,j) - u(i-1,j)) / (T[i]-T[i-1])*(t-T[i-1]);
		if (i>=n_dis)
			uu[j] = u(n_dis-1,j);
	}

	//cout << t << " " << 0 << " " << xx[0] << endl;
}


#ifdef TRANSWORHP_GRAPHICS
string ExplTransWorhp::setViewerTitle() {
    return string(" [ndis = " + std::to_string(n_dis) + ", multinodes = " + std::to_string(n_multinodes) + "]");
}

void ExplTransWorhp::updateViewer(DataStorage* ds, vector<double> &tmptime) {

	ds->setData(	twdiscretization->stuetzstellen(n_dis), tmptime.data(),
			block_x.data(), n_ode, n_ode,
			block_u.data(), n_ctrl,n_ctrl,
			G,  0,  0,
			lagrange_integral.data(), n_integral, n_integral);

	if (twparameter->showGrid) {
		ds->addDisData(n_dis, n_dis, T.data(), block_x.data(), block_u.data());
	}
}

void ExplTransWorhp::setTemptimeForViewer(vector<double> &tmptime) {

	tmptime = T;
}


#endif

int ExplTransWorhp::integrate_schrittweite(int index, int DGflag) {

#ifdef TW_WITH_BOOST
	if (twparameter->butchertableau >= 20 && twparameter->butchertableau < 30) {

		std::copy(pointer_x[index], pointer_x[index]+n_ode, tmp_ode_1.data());

		const int steps = butcher.integrate_with_boost(phase,T[index],T[index+1],tmp_ode_1,pointer_p, twparameter->butchertableau);

		for (int l = 0; l < n_ode; l++) {
			setx(index+1,l) = tmp_ode_1[l];
		}

		return steps;
	}
#endif

	if (butcher.stufen()==0) {

		double t  = T[index];
		double tend = T[index+1];

		std::copy(pointer_x[index],pointer_x[index]+n_ode,tmp_ode_1.data());
		std::copy(pointer_u[index],pointer_u[index]+n_ctrl,tmp_ctrl_1.data());
		std::copy(pointer_u[index+1],pointer_u[index+1]+n_ctrl,tmp_ctrl_2.data());

		int startflag = 1;

		const int steps = butcher.RungeKuttaF(phase,t,T[index],tend,tmp_ode_1.data(),tmp_ctrl_1.data(),tmp_ctrl_2.data(),pointer_p,startflag);

		startflag=2;

		for (int l=0; l<n_ode; l++) {
			setx(index+1,l) = tmp_ode_1[l];
		}

		return steps;
	}


	// Schrittweite
	double h  = butcher.h0;
	// Startpunkt (laeuft nach tend)
	double t  = T[index];
	// Endpunkt
	const double tend = T[index+1];

	//Anzahl der Integrationsschritte
	int steps = 0;

	/*
	for (int l=0; l<n_ode; l++)
		tmp_ode_1[l] = pointer_x[index][l];
	for (int l=0; l<n_ctrl; l++) {
		tmp_ctrl_1[l] = pointer_u[index][l];
		tmp_ctrl_2[l] = pointer_u[index+1][l];
	}
	*/
	std::copy(pointer_x[index],pointer_x[index]+n_ode,tmp_ode_1.data());
	std::copy(pointer_u[index],pointer_u[index]+n_ctrl,tmp_ctrl_1.data());
	std::copy(pointer_u[index+1],pointer_u[index+1]+n_ctrl,tmp_ctrl_2.data());

#if defined(TW_WITH_SUPERLU) || defined(WITH_LAPACK)
	// ROW-Verfahren (linear implizit)
	if (twparameter->butchertableau > 99) {

		steps = butcher->ROW(this, T[index], T[index+1], tmp_ode_1.data(), tmp_ctrl_1.data(), tmp_ctrl_2.data(), pointer_p);

		twcount_calls.integrate.f += steps;

		for (int l=0; l<n_ode; l++) {
			setx(index+1,l) = tmp_ode_1[l];
		}

		return steps;
	}
#endif

	if (DGflag == 0) {
		// speichert die benutzen Schrittweiten
		vector<double> auxSchrittweite;

		while (t<tend) {

			if (h+t>tend) {
				h = tend-t;
			}

			double auxH = h;

			int a = butcher.RungeKutta(phase,t,T[index],tend,tmp_ode_1.data(),tmp_ctrl_1.data(),tmp_ctrl_2.data(),pointer_p,h);

			steps++;

			if (a==-1) {
				cout << "Integrate: stepsize too small at " << t << ", stepsize:" << h << endl;
				break;
			}

			if (a==0) {
				if (t==tend) {
					for (int l=0; l<n_ode; l++) {
						setx(index+1,l) = tmp_ode_1[l];
					}
				}
				auxSchrittweite.push_back(auxH);
			}

			if (a==1) {
				//cerr << "erneut integrieren" << h << ", " << auxH << endl;
			}

		}

		schrittweiten.push_back(std::move(auxSchrittweite));

	} else {
		const int auxStepsize = butcher.getStepSize(); //Ordnung holen
		butcher.setStepSize(0); //damit das Verfahren als Verfahren ohne SW-Steuerung benutzt wird

		for (size_t m = 0; m < schrittweiten.at(index).size(); ++m) {

			h = schrittweiten[index][m];

			butcher.RungeKutta(phase,t,T[index],tend,tmp_ode_1.data(),tmp_ctrl_1.data(),tmp_ctrl_2.data(),pointer_p,h);

			steps++;
		}

		for (int l=0; l<n_ode; l++) {
			setx(index+1,l) = tmp_ode_1[l];
		}

		butcher.setStepSize(auxStepsize); // Ordnung setzen
	}


	return steps;
}

#ifdef _OPENMP
int ExplTransWorhp::integrate_schrittweite_parallel(int multi, int index, int DGflag) {

	if (butcher_parallel[multi].stufen()==0) {

		double t  = T[index];
		double tend = T[index+1];

		for (int l=0; l<n_ode; l++) {
			tmp_ode_parallel[index*n_ode+l] = pointer_x[index][l];
		}
		for (int l=0; l<n_ctrl; l++) {
			tmp_ctrl_1_parallel[index*n_ctrl+l] = pointer_u[index][l];
			tmp_ctrl_2_parallel[index*n_ctrl+l] = pointer_u[index+1][l];
		}

		int startflag = 1;

		int a = butcher_parallel[multi].RungeKuttaF(phase,t,T[index],tend,&tmp_ode_parallel[index*n_ode],&tmp_ctrl_1_parallel[index*n_ctrl],&tmp_ctrl_2_parallel[index*n_ctrl],pointer_p,startflag);
		startflag=2;

		for (int l=0; l<n_ode; l++) {
			setx(index+1,l) = tmp_ode_parallel[index*n_ode+l];
		}

		return a;
	}


	// Schrittweite
	double h  = butcher_parallel[multi].h0;
	// Startpunkt (laeuft nach tend)
	double t  = T[index];
	// Endpunkt
	const double tend = T[index+1];

	//Anzahl der Integrationsschritte
	int steps = 0;

	for (int l=0; l<n_ode; l++){
		tmp_ode_parallel[index*n_ode+l] = pointer_x[index][l];
	}
	for (int l=0; l<n_ctrl; l++) {
		tmp_ctrl_1_parallel[index*n_ctrl+l] = pointer_u[index][l];
		tmp_ctrl_2_parallel[index*n_ctrl+l] = pointer_u[index+1][l];
	}

#if defined(TW_WITH_SUPERLU) || defined(WITH_LAPACK)
	// ROW-Verfahren (linear implizit)
	if (twparameter->butchertableau > 99) {

		steps = butcher_parallel[multi]->ROW(this, T[index], T[index+1], &tmp_ode_parallel[index*n_ode], &tmp_ctrl_1_parallel[index*n_ctrl], &tmp_ctrl_2_parallel[index*n_ctrl], pointer_p);

		twcount_calls.integrate.f += steps;

		for (int l=0; l<n_ode; l++) {
			setx(index+1,l) = tmp_ode_parallel[index*n_ode+l];
		}

		return steps;
	}
#endif

	if (DGflag == 0) {
		// speichert die benutzen Schrittweiten
		vector<double> auxSchrittweite;

		while (t<tend) {

			if (h+t>tend) {
				h = tend-t;
			}

			double auxH = h;

			int a = butcher_parallel[multi].RungeKutta(phase,t,T[index],tend,&tmp_ode_parallel[index*n_ode],&tmp_ctrl_1_parallel[index*n_ctrl],&tmp_ctrl_2_parallel[index*n_ctrl],pointer_p,h);

			steps++;

			if (a==-1) {
				cout << "Integrate: stepsize too small at " << t << ", stepsize:" << h << endl;
				break;
			}

			if (a==0) {
				if (t==tend) {
					for (int l=0; l<n_ode; l++) {
						setx(index+1,l) = tmp_ode_parallel[index*n_ode+l];
					}
				}
				auxSchrittweite.push_back(auxH);
			}

			if (a==1) {
				//cerr << "erneut integrieren" << h << ", " << auxH << endl;
			}

		}

		schrittweiten[index] = auxSchrittweite;

	} else {

		const int auxStepsize = butcher_parallel[multi].getStepSize(); //Ordnung holen
		butcher_parallel[multi].setStepSize(0); //damit das Verfahren als Verfahren ohne SW-Steuerung benutzt wird

		for (size_t m = 0; m < schrittweiten.at(index).size(); ++m) {

			h = schrittweiten[index][m];

			butcher_parallel[multi].RungeKutta(phase,t,T[index],tend,&tmp_ode_parallel[index*n_ode],&tmp_ctrl_1_parallel[index*n_ctrl],&tmp_ctrl_2_parallel[index*n_ctrl],pointer_p,h);

			steps++;
		}

		for (int l=0; l<n_ode; l++) {
			setx(index+1,l) = tmp_ode_parallel[index*n_ode+l];
		}

		butcher_parallel[multi].setStepSize(auxStepsize); // Ordnung setzen
	}

	return steps;
}
#endif

void ExplTransWorhp::butcherInit(int index, double stepsize, bool verbose) {

	butcher.Init(index, stepsize, verbose);

#ifdef _OPENMP
	if (twparameter->parallel) {
		for (size_t i = 0; i < multinodes.size()-1; i++) {
			butcher_parallel[i].Init(index, stepsize, verbose);
		}
	}
#endif
}



void ExplTransWorhp::init0() {

	butcherInit(twparameter->butchertableau, twparameter->stepsize);

	n_var = n_ode*n_multinodes + n_ctrl * twdiscretization->stuetzstellen(n_dis) + n_param;
	n_con = n_ode*(n_multinodes-1);
	n_con += n_boxNeben*(twdiscretization->stuetzstellen(n_dis) - n_multinodes);

	if (n_con<0) n_con = 0;
	n_con += n_rand + n_neben * twdiscretization->stuetzstellen(n_dis);
	if (n_con<0) n_con = 0;

	tmp_gg_1 = vector<double>(n_con);
	tmp_gg_2 = vector<double>(n_con);


	phase->localinit();

	DS_obj.Init(1,         n_var);
	DS_ode.Init(n_ode,     n_ode+n_ctrl + n_param);
	DS_rand.Init(n_rand,   n_var);
	DS_neben.Init(n_neben, n_ode+n_ctrl + n_param);

	DS_integral.Init(n_integral, n_ode+n_ctrl + n_param);


	// Informationen ueber Optimierungsvariablen speichern

	int m = -1; // erste Punkt, muss Multiknoten sein!
	for (int i = 0; i < n_dis; ++i) {
		if (isMultinode(i)) {
			m++;
			for (int k = 0; k < n_ode; ++k) {
				typeX t;
				t.type = 0;
				t.n = k;
				t.multinodeL = m;
				t.dis = i;
				typeOfX.push_back(t);
			}
		}
		for (int k = 0; k < n_ctrl; ++k) {
			typeX t;
			t.type = 1;
			t.n = k;
			t.multinodeL = m;
			t.dis = i;
			typeOfX.push_back(t);
		}
	}
	for (int k = 0; k < n_param; ++k) {
		typeX t;
		t.type = 2;
		t.n = k;
		t.multinodeL = -1;
		t.dis = -1;
		typeOfX.push_back(t);
	}
}


void ExplTransWorhp::GetBoundaryIndices(vector<int> &indices, int d) {

	if ( X_low[x_index__(0,d)] == X_upp[x_index__(0,d)] ) {
		int ind = twdiscretization->stuetzstellen(1) ;
		indices.push_back(ind-1);
		//cout << it-twfolder->phases.begin() << " " << s<< " L " << ii << " " << ind-1 << endl;
	}
	//if ( X_low[x_index(n_dis-1,d)] == X_upp[x_index(n_dis-1,d)] ) {
	//	int ind = twdiscretization->stuetzstellen(n_dis) ;//+ ct;
	//	indices.push_back(ind-1);
	//cout << it-twfolder->phases.begin() << " " << s << " R " << ii << " " << ind-1 << endl;
	//}
}

}
