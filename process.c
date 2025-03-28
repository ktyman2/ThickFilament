// This file contains detailed processes.

void InitSingleProcessSubroutine(double *arr, int cnt) {
  double *ip;

  ip = arr;
  while (ip < arr + cnt) {
    *ip = 0.;
    ip++;
  }
}

// Initiate each single step 
void InitSingleProcess(void) {
  int n;

  currTimeStep++;
  // Forces are set to zero before SingleProcess().
  // forces of actins.
  InitSingleProcessSubroutine(act.f, (nActMe + nActCp) * NDIM);
  // Brownian forces of actins.
  InitSingleProcessSubroutine(act.fBr, (nActMe + nActCp) * NDIM);
  // forces of ACPs and motors.
  InitSingleProcessSubroutine(abp.f, (nAbpMe + nAbpCp) * NDIM);
  // Brownian forces of ACPs and motors.
  InitSingleProcessSubroutine(abp.fBr, (nAbpMe + nAbpCp) * NDIM);
  // forces actin on boundaries.
  if (bndMv.gTgl != 0 || rheoWay > 0) {
	for(n = 0; n < NDIM * 2 * NDIM; n++) { bnd.f[n] = 0.; }
  }
  if (actDgd.gTgl != 0) {
	actDgd.c = 0;
  }
  if (tglActFormDyn != 0) {
	cntNucAssAnn = 0;
  }

  if (currTimeStep == netForm.dur && motActiv.dur > 0) {
	Printf0("============================ Motor activation phase "
			"============================\n");
  }
  else if (currTimeStep == netForm.dur + motActiv.dur && rheo.dur > 0) {
	Printf0("=============================== Main measurement "
			"===============================\n");
  }
}

void SingleProcess(void) {
  // Initialize each step
  InitSingleProcess ();
  // Adjust size of subdomains depending on the number of particles in each 
  // Measure displacements of all particles to determine whether
  // neighboring list should be updated or not at a current time step.
  MeasureDisplacementForNL();
  // Calculate and apply stress or strain to a network for bulk 
  // rheology measurement
  if (rheoWay > 0 && currTimeStep >= netForm.dur + motActiv.dur) { 
	ApplyStressStrain(); 
  }
  // Repulsive forces between actin segments
  if (actF.rep.facStf > 0. || abpF.rep.facStf > 0.) {
	CalcRepulsiveForces(); 
  }
  // Spring forces for actin, ACP, and motor
  CalcSpringForces();
  // Bending forces for actin filaments
  if (actF.bend.facStf > 0.) {
	CalcFilaBendForces();
  }
  // Brownian forces for actin, ACP, and motor
  CalcBrownForces();
  // Two kinds of bending forces for ACPs
  if (!(motSA.gTgl != 0 && nAbp - nMot == 0)) {
	CalcAbpBendForces();
  }
  // Spring and bending forces acting on the backbone of multimerized motors
  if (motSA.gTgl != 0 && nMot > 0) { 
	CalcMotorBackboneForces();
  }

  /*--------------------------- Dynamics of actins ---------------------------*/
  // Polymerization of actin monomers on filaments
  if (actAss.tgl != 0 && actM.c > 0) {
	UpdateActinAssembly();
  }
  // Depolymerization of actin monomers on filaments
  if (actDis.tgl != 0 || actBst.tgl != 0) {
	UpdateActinDisassembly();
  }
  // Nucleation of actin
  if (actNuc.tgl != 0 && actM.c > 1) {
	UpdateActinNucleation();
  }
  // Severing of actin filaments
  if (actSev.tgl != 0) { 
	UpdateActinSevering();
  }
  // Eliminate expired elements in noActDyn.l
  if (tglActDyn != 0) {
	UpdateNoActinDynamicsList();
  }
  // If actin disassembly/severing is allowed with bound ABPs, this function 
  // should be executed here to let other adjacent subdomains know about ABP 
  // unbinding events caused by disassembly/severing.
  if (nCpu > 1 && tglActDisBstSev != 0) {
	UpdateActinAbpDynamicsEvents(0); 
  }
  /*----------------------- Dynamics of ACPs or motors -----------------------*/
  // Binding event of monomeric ACP or motor
  if ((acpMoBind.tgl != 0 && acpM.c > 0 && gTglImpAcpM != 0) 
		|| (motMoBind.tgl != 0 && motM.c > 0 && gTglImpMotM != 0)) {
	UpdateAbpMonomerBinding();
  }
  if (motSA.gTgl != 0) {
	if (motSA.to.gTgl != 0 && currTimeStep >= netForm.dur) {
		UpdateMotorTurnover();
	}
	if ((motM.c > 0 && gTglImpMotM != 0) || gTglImpMotM == 0) {
		UpdateMotorAssembly();
	}
  }
  // Unbinding event of inactive ACP or motor to monomer
  if (acpInaUnb.tgl != 0 || motInaUnb.tgl != 0) {
	UpdateInactAbpUnbinding();
  }
  // Binding event of inactive ACP or motor
  if (acpReb.tgl != 0 || motReb.tgl != 0) {
	UpdateAbpBinding();
  }
  // Unbinding event of active ACP or motor
  if (acpUnb.tgl != 0 || motUnb.tgl != 0) { 
	UpdateActiveAbpUnbinding();
  }
  // Walking event of active motor
  if (motWalk.tgl != 0) {
	UpdateMotorWalking();
  }
  // Eliminate expired elements in noAbpDyn.l 
  if (tglAbpAcInaDyn != 0 || acpMoBind.tgl != 0 || motMoBind.tgl != 0 
		|| tglActDisBstSevAbp != 0) { 
	UpdateNoAbpUnbindList();
  }
  // transfer actin severing/disassembly and ABP unbinding/binding events 
  // stored in sendAbpDyn.l to adjacent CPUs
  if (nCpu > 1 && tglAbpAcInaDyn != 0) { 
	UpdateActinAbpDynamicsEvents(1);   
  }
  if (updMono.tgl != 0) {
	if (nCpu > 1 && currTimeStep % updMono.prd  == 0  
			&& (tglAbpInaMoDyn != 0 || tglActMoDyn != 0)) {
		UpdateActinAbpMonomerList();
	}
  }
  if (motSA.gTgl != 0) {
	if (nCpu > 1 && currTimeStep % updMotN.prd == 0 && updMotN.tgl != 0) {
		UpdateMotorNucleCounter();
	}
  }
  /*-------------------------------- Boundary --------------------------------*/
  // Unbinding and binding of filaments from boundaries
  if (bnd.gotF != 0) {
	if (bndUnb.gTgl != 0) {
		UpdateBoundaryActinUnbindMature();
	}
	if (bndReb.gTgl != 0) {
		UpdateBoundaryActinBinding();
	}
  }
  /*--------------------------------------------------------------------------*/
  // Calculate new positions of all particles based on forces 
  UpdateNewLocation();
  // Apply the periodic boundary condition to all particles
  ApplyBoundCondAll();
  // Only for a case with more than 1 CPU
  if (nCpu > 1) { 
	// Adjust the locations of boundaries following stress-strain relationship
	if (bndMv.gTgl != 0 && currTimeStep >= netForm.dur 
			+ pres.dur + motActiv.dur) {
		UpdateBoundaryLocation();		
	}
	// Adjust size of subdomains depending on the number of particles in each 
	// subdomain. The aim is to maintain the number at similar level.
	// MoveParticles() and CopyParticles() should be excuted right after 
	// UpdateSubdomSectionLocation().
	if (updSubdSize.tgl != 0) {
		if ((currTimeStep == netForm.dur && gTglLoadNetData == 0) 
				|| (currTimeStep % updSubdSize.prd == 0
				&& !(currTimeStep < netForm.dur && gTglLoadNetData == 0))) {
			UpdateSubdomSectionLocation(0);
			if (rank == 0) { 
				RecordSubdomSectionLocation(); 
			}
			MoveParticles();
			CopyParticles();
			UpdateNeighborList();
		}
	}
	UpdateLongChainNormal(); 
	// Move particles escaping from a current subdomain to adjacent subdomains
	MoveParticles();
	// Copy particles located near boundaries to adjacent subdomains 
	// for synchronization
	CopyParticles();
  }
}


void InitToggleParameters(int mode) {
  // Actin assembly or disassembly is on or off 
  // depending on parameters in 'condition'.
  if (nAct > 0 && (mode == 2 || (mode == 1 && gTglActDynPres != 0)
		|| (mode == 0 && gTglActDynNF != 0))) {
	actAss.tgl = actAss.gTgl;
	actDis.tgl = actDis.gTgl;
	actBst.tgl = actBst.gTgl;
	actSev.tgl = actSev.gTgl;
	actAnn.tgl = actAnn.gTgl;
	actNuc.tgl = actNuc.gTgl;
	actCap.tgl = actCap.gTgl;
	actUnc.tgl = actUnc.gTgl;
  }
  else {
	actAss.tgl = 0;
	actDis.tgl = 0;
	actBst.tgl = 0;
	actSev.tgl = 0;
	actAnn.tgl = 0;
	actNuc.tgl = 0;
	actCap.tgl = 0;
	actUnc.tgl = 0;
  }
  // Binding of monomeric ACPs or motors is on or off.
  if ((mode == 2 && nAbp > 0) || (mode == 1 && ((gTglAcpDynPres != 0 
		&& nAbp - nMot > 0) || (gTglMotUnbRebPres != 0 && nMot > 0)))
		|| (mode == 0 && ((gTglAcpDynNF != 0 && nAbp - nMot > 0) 
		|| (gTglMotUnbRebNF != 0 && nMot > 0)))) { 
	acpMoBind.tgl = acpMoBind.gTgl;
	motMoBind.tgl = motMoBind.gTgl;
  }
  else { 
	acpMoBind.tgl = 0; 
	motMoBind.tgl = 0; 
  }
  // Unbinding of active ACPs or motors and 
  // unbinding and binding of inactive ACPs of motors are on or off.
  if (nAbp - nMot > 0 && (mode == 2 || (mode == 1 && gTglAcpDynPres != 0)
		|| (mode == 0 && gTglAcpDynNF != 0))) {
	acpUnb.tgl = acpUnb.gTgl;
	acpReb.tgl = acpReb.gTgl;
	acpInaUnb.tgl = acpInaUnb.gTgl;
  }
  else {
	acpUnb.tgl = 0;
	acpReb.tgl = 0;
	acpInaUnb.tgl = 0;
  }
  // Unbinding of active motors and the unbinding or binding of inactive
  // motors are on or off. 
  if (nMot > 0 && (mode == 2 || (mode == 1 && gTglMotUnbRebPres != 0)
		|| (mode == 0 && gTglMotUnbRebNF != 0))) {
	motUnb.tgl = motUnb.gTgl;
	motReb.tgl = motReb.gTgl;
	motInaUnb.tgl = motInaUnb.gTgl;
  }
  else {
	motUnb.tgl = 0;
	motReb.tgl = 0;
	motInaUnb.tgl = 0;
  }
  // Walking of active motors is on or off.
  if (nMot > 0 && (mode == 2 || (mode == 1 && gTglMotWalkPres != 0)
		|| (mode == 0 && gTglMotWalkNF != 0))) {
	motWalk.tgl = motWalk.gTgl;
  }
  else { 
	motWalk.tgl = 0; 
  }
  // Toggle parameters representing multiple toggles are on or off.
  tglAcpAcInaDyn = (acpUnb.tgl != 0 || acpReb.tgl != 0 
		|| acpInaUnb.tgl != 0) ? 1 : 0;
  tglMotAcInaDyn = (motUnb.tgl != 0 || motReb.tgl != 0 || motWalk.tgl != 0 
		|| motInaUnb.tgl != 0) ? 1 : 0;
  tglAbpAcInaDyn = (tglAcpAcInaDyn != 0 || tglMotAcInaDyn != 0) ? 1 : 0;
  tglAbpInaMoDyn = (acpMoBind.tgl != 0 || motMoBind.tgl != 0
		|| acpInaUnb.tgl != 0 || motInaUnb.tgl != 0) ? 1 : 0;
  tglActDisBstSev = (actBst.tgl != 0 || actDis.tgl != 0 
		|| actSev.tgl != 0) ? 1 : 0;
  tglActDisBstSevAbp = ((actDis.tgl != 0 && actDis.facKWA > 0.) 
		|| (actBst.tgl != 0 && actBst.facKWA > 0.) 
		|| (actSev.tgl != 0 && actSev.facKWA > 0.)) ? 1 : 0;
  tglActFormDyn = (actAss.tgl != 0 || actNuc.tgl != 0 || actAnn.tgl != 0) 
		? 1 : 0;
  tglActMoDyn = (actAss.tgl != 0 || actDis.tgl != 0 || actBst.tgl != 0 
		|| actNuc.tgl != 0) ? 1 : 0;
  tglActDyn = (tglActMoDyn != 0 || actSev.tgl != 0 || actAnn.tgl != 0) ? 1 : 0;
  tglMbDyn = 0;
  tglMbNucDyn = 0;
  tglMbNucLocAre = 0;
  tglNeiAbpSC = (acpReb.gTgl != 0 || motReb.gTgl != 0 
		|| abpF.rep.facStf > 0.) ? 1 : 0;
  tglNeiAbpDC = (abpF.rep.facStf > 0.) ? 1 : 0;
}

// Main process
// There are two separate parts for each rheology measurement.
// rheoWay = 0 -> segment-rheology, 1 -> bulk rheology
void MainProcess(void) {
  // Synchronization
  MPI_Barrier(MPI_COMM_WORLD);
  // Initial output.
  if (rank == 0) { 
	InitOutput(); 
  }
  InitToggleParameters(0);
  // This loop is iterated until currTimeStep reaches a last time step.
  // "netForm.dur" is the time step required for network formation.
  // "pres.dur" is the time step required to apply prestrain/prestress.
  // "rheo.dur" is the time of rheological measurement in "condition".
  while (currTimeStep <= netForm.dur + rheo.dur + pres.dur + motActiv.dur) {
	if (currTimeStep == netForm.dur) {
		InitToggleParameters(1);
	}
	if (currTimeStep == netForm.dur + pres.dur + motActiv.dur) {
		InitToggleParameters(2);
	}
	// Excute the main part of the process run every time step.
	SingleProcess();
	// Check and record things after SingleProcess().
	TaskAfterSingleProcess();
  }
  Printf0("All process is done!\n"); 
  // Finalize MPI processes
  MPI_Finalize();
}

void TaskAfterSingleProcess(void) {
  char fn[80];

  /*------------------------------ Record progress ---------------------------*/
  // Record progress information and the mechanical energy of networks.
  if (currTimeStep % recProg.prd == 0 || currTimeStep == 1) {
    RecordProgress(); 
  }
  /*------------------------- Rheological measurement ------------------------*/
  if (recSecStre.tgl != 0) {
	// Measure internal elastic and viscous stress
	RecordElasViscStress(recSecStre.prd, nSecStreDiv, recSecStre.mode);
  }
  /*----------------------- Boundaries and subdomains ------------------------*/
  if ((bndUnb.gTgl != 0 || bndReb.gTgl != 0) && recBndUnbReb.tgl != 0) {
	if (currTimeStep % recBndUnbReb.prd == 0) {
		RecordBoundaryActinUnbindBind();
	}
  }
  if (bndMv.gTgl != 0 && recBndLoc.tgl != 0 && rank == 0) {
	if (currTimeStep % recBndLoc.prd == 0) {
		RecordBoundaryLocation();
	}
  }
  /*--------------------------- Network information --------------------------*/
  if (recConf.tgl != 0) {
	// Record network information including the chains and positions of 
	// particles.
	if (currTimeStep % recConf.prd == 0 || currTimeStep == 1) {
		RecordConfig(GenFileName("AccuConf"), recConf.prd);
	}
  }
  if (recConfVmd.tgl != 0) {
	// Record networks for visualization via VMD. 
	if (currTimeStep % recConfVmd.prd == 0 || currTimeStep == 1) {
		RecordConfigVmd(0);
	}
  }
  /*--------------------------- Energy calculation ---------------------------*/
  if (recE.tgl != 0) {
	if (currTimeStep % recE.prd == 0 || currTimeStep == 1) {
	    RecordMechEnergy(0);
		if (recPerc.tgl != 0) { 
			RecordMechEnergy(1); 
		}	
		if (findSupp.tgl != 0) {
			RecordMechEnergy(2);
		}
	}
  }
  /*-------- Information in unit of elements, segments, and filaments --------*/

  if (recInfo.tgl != 0) { 
	if (currTimeStep % recInfo.prd == 0) {
		RecordActinInfo(recInfo.prd);
		RecordIndvSegFilaInformation(recInfo.prd);
	}
  }
  if (confVmdInfo.tgl != 0) {
	// If a specific coloring method is used, the accumulated values of chain  
	// lengths and forces have to recorded.
    RecordAccuLengthForces();
	if (currTimeStep % confVmdInfo.prd == 0) {
		ResetAccuLengthForces();
	}
  }
  if (recLongF.tgl != 0) {
	RecordAccuLongSpringForces(recLongF.prd);
  }
  /*------------------ Dynamics behaviors of actin and ABPs ------------------*/

  /*-------------------- Structural properties of network --------------------*/
  if (recCrsDist.tgl != 0) {
	if (currTimeStep % recCrsDist.prd == 0 || currTimeStep == 1) {
		RecordCrosslinkDistance();
	}
  }
  if (recFilaL.tgl != 0) {
	if (currTimeStep % recFilaL.prd == 0) {
		RecordFilamentLength(recFilaL.prd);
	}
  }

  // If an error signal appears, the run shoulde be terminated.
  if (stopSig != 0) {
	printf("Due to the large force, the code is stopped at currTimeStep = "
		"%lld!\n", currTimeStep);
	while(1);
		exit(-1);
	}
}

void PrepareStateWoNetworkData(void) {
  int n, nIFilaP;

  nAct = 0;
  nActMe = 0;
  nActCp = 0;
  actM.c = 0; 
  nActFilaMe = 0;

  nAbp = 0;
  nAbpMe = 0;
  nAbpCp = 0; 
  if (gTglImpAcpM != 0) { acpM.c = 0; }
  if (gTglImpMotM != 0) { motM.c = 0; }

  recTraj.actMe.c = 0; 
  recTraj.abpMe.c = 0; 
  if (recTraj.gTglCho == 0 && (rheoWay == 0 || rheoWay == 2)) {
	Printf0("Warning: new components won't be chosen for tracking trajectory "
			"without loaded network data. Nothing will be traced!\n\n");
  }

  if (rheoWay > 0) { 
	meaStreParMe.c = 0;
	appStraParMe.c = 0;
  }
  neigh.c = 0;
  if (actSev.gTgl != 0 || actNuc.gTgl != 0 || actAnn.gTgl != 0) {
	nIFilaP = (int)(nActGoal / nCpu);
	for(n = 0; n < nIFilaP; n++) {
		iFilaP.l[n] = rank * nIFilaP + n;
	}
	iFilaP.c = nIFilaP;
  }
  UpdateChainList();
  
  if (rheoWay > 0) {
	if (bndReb.gTgl == 0) {
		Printf0("Error: if Config is not loaded, the binding of actins "
				"on boundaries should be allowed for bulk rheology!\n\n");
		exit(-1);
	}
  }
}

void PrepareNotBulkRheology(void) {
  int n, k;
  double *pR;
  char dirCh[] = "xyz";

  if (rheoWay == 0) {
	if (recTraj.gTglCho == 0 && recTraj.act.c == 0 && recTraj.abp.c == 0) {
		Printf0("Warning: there is no information of actins or ABPs whose "
				"trajectories are recorded in Config. Nothing will be "
				"traced!\n\n");
	}
  }
  FOR_NDIM(k) {
	if (confPbc[k] != 0) {
		if (pbc[k] == 0) {
			for(n = 0; n < nAct + nAbp; n++) {
				pR = (n < nAct) ? &P2(rAct,n,k) : &P2(rAbp,n - nAct,k);
				if (*pR >= rGrid[k][nGrid[k] - 1]) { *pR -= dimDom[k]; }
				else if (*pR < rGrid[k][0]) { *pR += dimDom[k]; }
			}
			SeverFilaOnPlane(k, 1);
		}
	}
	else {
		Printf0("Warning: the loaded network doesn't have a periodic boundary "
				"condition in the %c direction!\n\n", dirCh[k]);
	}
  }
}
