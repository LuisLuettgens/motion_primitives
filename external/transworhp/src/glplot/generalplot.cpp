#include "generalplot.h"

#include "xopt_eps.h"

namespace tw {

GeneralPlot::GeneralPlot(char c_, int d, int ind) : BasePlot(ind), data(c_,d) {

}

GeneralPlot::~GeneralPlot() {}

int GeneralPlot::GetDgl() const {
	return data.n;
}

int GeneralPlot::rawHighLow(DataStorage &ds) {

	for (int i = 0; i < ds.getLength(); i++) {
		double v = ds.getData(data,i) * scaledata;
		if (v < Low || i==0) {
			Low = v;
		}
		if (v > High || i==0) {
			High = v;
		}
	}

	return 0;
}


void GeneralPlot::SetMaxTime(DataStorage::TimeMode_e timemode, double time) {
	High2 = time;
	if (timemode == DataStorage::index_as_time) {
		High2 = time-1;
	}
}

void GeneralPlot::SetMinTime(DataStorage::TimeMode_e timemode, double time) {
	Low2 = time;
	if (timemode == DataStorage::index_as_time) {
		Low2 = time-1;
	}
}

void GeneralPlot::drawData(DataStorage &ds, DataStorage &/*dstop*/, int ii) {

	// Speicher fuer Plot
	std::vector<GLfloat> plotData(ds.getLength()*2);

	// Speicher fuer fixe Randwerte
	std::vector<GLfloat> plotDataPoints;
	size_t plotDataPointsSize = 0;
	if (dot_indices.size() > static_cast<size_t>(ii)) {
		plotDataPointsSize = dot_indices[ii].size();
	}
	plotDataPoints.resize(plotDataPointsSize*2);

	// wenn diff = 1: Gauss-pmTW Verfahren
	int diff = 0;
	if (ds.getDisLengthX() != ds.getDisLengthU()) {
		diff = std::abs(ds.getDisLengthX() - ds.getDisLengthU())-1;
	}

	// Speicher fuer Zeit-Diskretisierung
	std::vector<GLfloat> plotTimeDis;
	size_t plotTimeDisSize = ds.getDisLengthX()-diff*2;
	plotTimeDis.resize(plotTimeDisSize*4);

	// Speicher fuer Zeit-Diskretisierung-Punkte
	std::vector<GLfloat> plotTimeDisPoints;
	plotTimeDisPoints.resize(plotTimeDisSize*2);

	glLineWidth(1.5f);
	DrawCompareCurve(ds);
	DrawDynCompareCurve(ds);

	switch (data.c) {
		case 'x':
			SetColor(TWcolor::Black);
			break;
		case 'u':
			SetColor(TWcolor::Red);
			break;
		case 'l':
			SetColor(TWcolor::Cyan);
			break;
	}

	glPushMatrix();

	// Trafo: ehemals x: BasePlot::MapTimeToX(..); y: BasePlot::MapToY(..);
	if (!IsIcon()) {
		glTranslatef(LBorder, TBorder, 0.0f);
		if (ds.timemode == DataStorage::time_is_time) { //feste Endzeit
			glScalef((GetWidth()-LBorder-RBorder)/(ds.MaxTF-ds.MinT0), scaledata*yscale, 1.0f);
			glTranslatef(-ds.MinT0, -Low, 0.0f);
		} else if (ds.timemode == DataStorage::fixed_float_as_time) { // freie Endzeit
			glScalef((GetWidth()-LBorder-RBorder), scaledata*yscale, 1.0f);
			glTranslatef(ds.start, -Low, 0.0f);
			glScalef(ds.getTF()/ds.MaxTF, 1.0f, 1.0f);
		}
	} else {
		if (ds.timemode == DataStorage::time_is_time) { //feste Endzeit
			glScalef(GetWidth()/(ds.MaxTF-ds.MinT0), scaledata*yscale, 1.0f);
			glTranslatef(-ds.MinT0, -Low, 0.0f);
		} else if (ds.timemode == DataStorage::fixed_float_as_time) { // freie Endzeit
			glScalef((GetWidth()), scaledata*yscale, 1.0f);
			glTranslatef(ds.start, -Low, 0.0f);
			glScalef(ds.getTF()/ds.MaxTF, 1.0f, 1.0f);
		}
	}

	if (dot_indices.size() > static_cast<size_t>(ii)) {
		int i = 0;
		for (auto it = dot_indices[ii].cbegin(); it != dot_indices[ii].end(); ++it) {
			plotDataPoints[i*2] = static_cast<GLfloat>( ds.getData(Selector('t',0),*it) );
			plotDataPoints[i*2+1] = static_cast<GLfloat>( ds.getData(data,*it) );
			i++;
		}
	}

	for (int i = diff, j = 0; i < ds.getDisLengthX()-diff; i++, j++) {

		const GLfloat aux = static_cast<GLfloat>(ds.getData(Selector('T',0),i));

		plotTimeDis[j*4] = static_cast<GLfloat>( aux );
		plotTimeDis[j*4+1] = static_cast<GLfloat>( Low );
		plotTimeDis[j*4+2] = static_cast<GLfloat>( aux );
		plotTimeDis[j*4+3] = static_cast<GLfloat>( High );
	}

	if (data.c == 'x') {
		for (int i = diff, j = 0; i < ds.getDisLengthX()-diff; i++,j++) {
			plotTimeDisPoints[j*2] = static_cast<GLfloat>(ds.getData(Selector('T',0),i));
			plotTimeDisPoints[j*2+1] = static_cast<GLfloat>(ds.getData(Selector('X',data.n),i-diff));
		}
	} else if (data.c == 'u') {
		for (int i = 0; i < ds.getDisLengthU(); i++) {
			plotTimeDisPoints[i*2] = static_cast<GLfloat>(ds.getData(Selector('T',0),i+diff));
			plotTimeDisPoints[i*2+1] = static_cast<GLfloat>(ds.getData(Selector('U',data.n),i));
		}
	}

	for (int i = 0; i < ds.getLength(); i++) {
		plotData[i*2] = static_cast<GLfloat>( ds.getData(Selector('t',0),i) );
		plotData[i*2+1] = static_cast<GLfloat>( ds.getData(data,i) );
	}

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glEnableClientState(GL_VERTEX_ARRAY);

	glLineWidth(1.5f);

	glVertexPointer(2, GL_FLOAT, 0, plotData.data());
	glDrawArrays(GL_LINE_STRIP, 0, ds.getLength());

	glLineWidth(1.0f);

	if (dot_indices.size() > static_cast<size_t>(ii)) {
		glPointSize(5.0f);
		glVertexPointer(2, GL_FLOAT, 0, plotDataPoints.data());
		glDrawArrays(GL_POINTS, 0, plotDataPointsSize);
		glPointSize(1.0f);
	}

	// Zeit-Diskretisierung
	if (ds.getDisLengthX()) {

		glPointSize(5.0f);
		glVertexPointer(2, GL_FLOAT, 0, plotTimeDisPoints.data());
		glDrawArrays(GL_POINTS, 0, ds.getDisLengthX()-diff*2);
		glPointSize(1.0f);

		SetColor(TWcolor::Cyan);
		glVertexPointer(2, GL_FLOAT, 0, plotTimeDis.data());
		glDrawArrays(GL_LINES, 0, (ds.getDisLengthX()-diff*2)*2);
	}

	glDisableClientState(GL_VERTEX_ARRAY);

	glPopMatrix();
}


bool GeneralPlot::MouseInput(DataStorage &ds, int /*button*/, const Point<int> &p) {

	if (hasControlData()) {

		const double xx = MapFromX(ds, p.x-GetPos().x);
		const int ctrlx = static_cast<int>(ds.getLength()*xx);

		const double ctrly = MapFromY(p.y-GetPos().y);
		ds.setData(data,ctrlx,ctrly);

		return true;
	}

	return false;
}


void GeneralPlot::epsData(DataStorage &ds, DataStorage &dstop, EpsWriter *epsw) {

	epsw->SetLineColor(TWcolor::Red);// "col2";

	for (int i=0;i<ds.getLength();i++) {
		epsw->AddPoint(epsw->epsMapTimeToX(ds,i),
		               epsw->epsMapToY(ds.getData(data,i)*scaledata));
	}
	if (!epsw->absetzen)
		epsw->Line();
	else {
		epsw->SetLineColor(TWcolor::Green);
		epsw->Dot(static_cast<float>(ds.getLength()));
		epsw->LineWidth(0.6f);
	}
	epsw->ClearPoints();


	if (dstop.getLength()) {
		epsw->SetLineColor(TWcolor::Green);// "col3";
		epsw->LineWidth(1.2f);

		for (int i=0;i<dstop.getLength();i++) {
			epsw->AddPoint(epsw->epsMapStopTimeToX(ds,dstop,i),
			               epsw->epsMapToY(dstop.getData(data,i)*scaledata));
		}
		epsw->Line();
		epsw->ClearPoints();

	}
}


double GeneralPlot::MapFromX(DataStorage &ds, short x) const {

	int width = GetWidth();

	if (ds.timemode==DataStorage::index_as_time) {
		float ww = width-LBorder-RBorder;
		float w2 = LBorder;
		int a = ds.getTF();
		return ((x-w2)/ww*(a-1.)+.5)/a;
	} else {
		return BasePlot::MapFromX(x);
	}
}

}
