/*
 * ntupleQCD.h
 *
 *  Created on: 24 Aug 2016
 *      Author: jkiesele
 */

#ifndef ntupleQCD_H_
#define ntupleQCD_H_

#include "interface/basicAnalyzer.h"
#include "interface/sampleCollection.h"
#include "classes/DelphesClasses.h"


class ntupleQCD: public d_ana::basicAnalyzer{
public:
	ntupleQCD():d_ana::basicAnalyzer(){}
	~ntupleQCD(){}


private:
	void analyze(size_t id);

	void postProcess();
};





#endif /* ntupleQCD_H_ */
