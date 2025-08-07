#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_BCRec.H>
//#include <AMReX_Interpolater.H>

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MLABecLaplacian.H>

//#include <mechanism.H>
#include <PelePhysics.H>

using namespace amrex;

pele::physics::PeleParams<pele::physics::transport::TransParm<
  pele::physics::PhysicsType::eos_type,
  pele::physics::PhysicsType::transport_type>>
  trans_parms;

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile infile=f1 [options] \n\tOptions:\n";
  exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    if (argc < 2)
      print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc,argv);

    if (pp.contains("verbose"))
      AmrData::SetVerbose(true);

    std::string plotFileName; pp.get("infile",plotFileName);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
      // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();

    trans_parms.initialize();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    Vector<std::string> spec_names;
    pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);
    std::map<std::string, int> nametoindex;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      nametoindex[spec_names[n]] = n;
    }    
    std::string fuelName = "H2"; pp.query("fuelName",fuelName);
    const int fuelIndex = nametoindex[fuelName];   
    std::string progressVar = "Y("+fuelName+")";
    ParmParse pptrans("transport");
    int do_soret = 0; pptrans.query("use_soret",do_soret);
    int do_wbar = 1; pp.query("use_wbar",do_wbar);
    int idYin = -1;
    int idTin = -1;
    int idRin = -1;
    int idCin = -1;
    
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    const std::string spName= "Y(" + spec_names[0] + ")";
    const std::string TName = "temp";
    const std::string RName = "density";
    const std::string CName = "I_R(" + fuelName + ")";
    
    for (int i=0; i<plotVarNames.size(); ++i)
    {
      if (plotVarNames[i] == spName) idYin = i;
      if (plotVarNames[i] == TName)  idTin = i;
      if (plotVarNames[i] == RName)  idRin = i;
      if (plotVarNames[i] == CName)  idCin = i;
    }
    if (idYin<0 || idTin<0 || idRin<0 || idCin<0)
      Abort("Cannot find required data in pltfile");

    const int nCompIn  = NUM_SPECIES+4; //all species, temperature, density, chemical source term for fuel
    const int nCompOut = 6; //each component of sd, plus sd, plus sdtilde
    Vector<std::string> outNames(nCompOut);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    const int idYlocal = 0; // Ys start here
    const int idClocal = NUM_SPECIES;
    const int idTlocal = NUM_SPECIES+1;   // T starts here
    const int idRlocal = NUM_SPECIES+2; // R starts here
    for (int i=0; i<NUM_SPECIES; ++i) {
      destFillComps[i] = idYlocal + i;
      inNames[i] =  "Y(" + spec_names[i] + ")";
    }
    destFillComps[idClocal] = idClocal;
    destFillComps[idTlocal] = idTlocal;
    destFillComps[idRlocal] = idRlocal;
    inNames[idTlocal] = TName;
    inNames[idRlocal] = RName;
    inNames[idClocal] = CName;
    outNames[0] = "Sd("+fuelName+")_dY";
    outNames[1] = "Sd("+fuelName+")_dT";
    outNames[2] = "Sd("+fuelName+")_dW";
    outNames[3] = "Sd("+fuelName+")_C";
    outNames[4] = "Sd("+fuelName+")";
    outNames[5] = "Sd("+fuelName+")tilde";

    Real rhou; pp.get("rhou",rhou);
    int idProglocal=-1;
    for (int i = 0; i < nCompIn; i++) {
      if (progressVar == inNames[i]) {
	idProglocal=i;
	break;
      }
    }
    if (idProglocal == -1) {
      Abort("Progress variable not being loaded");	    
    }

    Vector<int> sym_dir(AMREX_SPACEDIM,0);
    pp.queryarr("sym_dir",sym_dir,0,AMREX_SPACEDIM);

    Vector<int> is_per(AMREX_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Print() << is_per[idim] << " ";
    }
    Print() << "\n";
    BCRec gradVarBC;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        gradVarBC.setLo(idim,BCType::foextrap);
        gradVarBC.setHi(idim,BCType::foextrap);
        if ( is_per[idim] ) {
            gradVarBC.setLo(idim, BCType::int_dir);
            gradVarBC.setHi(idim, BCType::int_dir);
        }
    }

    int coord = 0;

    amrex::RealBox real_box({AMREX_D_DECL(amrData.ProbLo()[0], amrData.ProbLo()[1], amrData.ProbLo()[2])},
                            {AMREX_D_DECL(amrData.ProbHi()[0], amrData.ProbHi()[1], amrData.ProbHi()[2])});

    Vector<MultiFab> coeffMF(Nlev); //holds transport coefficients
    Vector<MultiFab> gradMF(Nlev); //holds all the gradients, need gradY, gradW and gradT
    Vector<MultiFab> grad2MF(Nlev); //holds all the second derivatives
    Vector<MultiFab> deriveMF(Nlev); //just holds mean molecular weight   
    Vector<Geometry> geoms(Nlev);
    Vector<BoxArray> grids(Nlev);
    Vector<DistributionMapping> dmap(Nlev);
    Vector<MultiFab> outdata(Nlev);
    Vector<MultiFab> indata(Nlev);
    const int nGrow = 1;
    const int nCoeffs = 3; //standard, soret and molecular weight coeffs for H2
    const int nGrads = 3*AMREX_SPACEDIM; //gradFuel gradWbar gradT
    const int n2Grads = 3*AMREX_SPACEDIM*AMREX_SPACEDIM; //each diffusive flux has D dims
    const int nDerive = 1 + 3*AMREX_SPACEDIM; //Wbar and the 3 different diffusive fluxes 
    // Read data on all the levels                                                                                                   
    for (int lev=0; lev<Nlev; ++lev) {
      const BoxArray ba = amrData.boxArray(lev);
      grids[lev] = ba;
      dmap[lev] = DistributionMapping(ba);
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&real_box,coord,&(is_per[0]));
      indata[lev].define(grids[lev], dmap[lev], nCompIn, nGrow);
      coeffMF[lev].define(grids[lev], dmap[lev], nCoeffs, nGrow);
      deriveMF[lev].define(grids[lev], dmap[lev], nDerive, nGrow);
      gradMF[lev].define(grids[lev],dmap[lev], nGrads, nGrow);
      grad2MF[lev].define(grids[lev],dmap[lev], n2Grads, nGrow);
      outdata[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      
      Print() << "Reading data for level: " << lev << std::endl;
      amrData.FillVar(indata[lev], lev, inNames, destFillComps);
      indata[lev].FillBoundary(0,1,geoms[lev].periodicity());
    }

    for (int lev=0; lev<Nlev; ++lev)
    {
      //First compute transport properties and specific enthalpy
      // Get the transport data pointer
      auto const* ltransparm = trans_parms.device_parm();
      
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (amrex::MFIter mfi(indata[lev], amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {

        const Box& bx = mfi.tilebox();
        Array4<Real const> const& Y_a = indata[lev].const_array(mfi,idYlocal);
        Array4<Real const> const& T_a = indata[lev].const_array(mfi,idTlocal);
        Array4<Real const> const& rho_a = indata[lev].const_array(mfi,idRlocal);
	
	Array4<Real> const& w_a = deriveMF[lev].array(mfi,0);
	Array4<Real> const& spec_coeff_a = coeffMF[lev].array(mfi,0);
	Array4<Real> const& soret_coeff_a = coeffMF[lev].array(mfi,1);
	Array4<Real> const& molar_coeff_a = coeffMF[lev].array(mfi,2);
	amrex::ParallelFor(bx, [=]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
	  
          Real Yloc[NUM_SPECIES] = {0.0};
	  Real Dcoeff[NUM_SPECIES] = {0.0};
	  Real chi[NUM_SPECIES] = {0.0};
	  Real imw[NUM_SPECIES] = {0.0};
	  for (int n=0; n<NUM_SPECIES; ++n) {
	    Yloc[n] = Y_a(i,j,k,n);
          }

	  //Dummies unused due to falses
	  Real lam,mu_dummy,xi_dummy,mmw;
	  //get D, lam, and chi
	  const bool get_xi = false;
	  const bool get_mu = false;
	  const bool get_lam = false;
	  const bool get_Ddiag = true;
	  const bool get_chi = (do_soret == 1);
	  pele::physics::transport::SimpleTransport::transport(get_xi,get_mu,get_lam,get_Ddiag,get_chi,T_a(i,j,k),rho_a(i,j,k),Yloc,Dcoeff,chi,mu_dummy,xi_dummy, lam,ltransparm);

	  //inverse molecular weights
	  get_imw(imw);
	  //mean molecular weight
	  CKMMWY(Yloc,mmw);
	  mmw *= 1.0e-3; //cgs->mks
	  w_a(i,j,k) = mmw;
	  for (int n=0; n< NUM_SPECIES; n++) {
	    imw[n] *= 1.0e3; //cgs->mks
	    Dcoeff[n] *= 1.0e-1; //cgs->mks
	    chi[n] *= Dcoeff[n]; //chi->theta
	  }
	  
	  spec_coeff_a(i,j,k) = Dcoeff[fuelIndex] * mmw * imw[fuelIndex];
	  molar_coeff_a(i,j,k) = (do_wbar == 1) ? Dcoeff[fuelIndex] * Yloc[fuelIndex] * imw[fuelIndex] : 0.0;
	  soret_coeff_a(i,j,k) = (do_soret == 1 && fuelName == "H2") ? 0.664*chi[fuelIndex]/T_a(i,j,k) : 0.0; 

	});
      }

      Print() << "Coefficients, enthalpy and MMW derived for level " << lev << std::endl;
    }


    // Get face-centered gradients from MLMG                                                                                          
    LPInfo info;
    info.setAgglomeration(1);
    info.setConsolidation(1);
    info.setMetricTerm(false);
    info.setMaxCoarseningLevel(0);
    MLPoisson poisson({geoms}, {grids}, {dmap}, info);
    poisson.setMaxOrder(4);
    std::array<LinOpBCType, AMREX_SPACEDIM> lo_bc;
    std::array<LinOpBCType, AMREX_SPACEDIM> hi_bc;
    for (int idim = 0; idim< AMREX_SPACEDIM; idim++){
       if (is_per[idim] == 1) {
          lo_bc[idim] = hi_bc[idim] = LinOpBCType::Periodic;
       } else {
          if (sym_dir[idim] == 1) {
             lo_bc[idim] = hi_bc[idim] = LinOpBCType::reflect_odd;
          } else {
             lo_bc[idim] = hi_bc[idim] = LinOpBCType::Neumann;
          }
       }
    }
    poisson.setDomainBC(lo_bc, hi_bc);


    // Need to apply the operator to ensure CF consistency with composite solve                                                       
    int nGrowGrad = 0;                   // No need for ghost face on gradient   

    Vector<Array<MultiFab,AMREX_SPACEDIM> > grad(Nlev);
    Vector<std::unique_ptr<MultiFab>> phi;
    Vector<MultiFab> laps;
    
    //for (int n =0; n<NUM_SPECIES; n++) {
    Print() << "Calculating grad"+inNames[idProglocal] << std::endl;
    for (int lev = 0; lev < Nlev; ++lev) {
      for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
	const auto& ba = grids[lev];
	grad[lev][idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
			       dmap[lev], 1, nGrowGrad);
      }
      phi.push_back(std::make_unique<MultiFab>(indata[lev],amrex::make_alias,idProglocal,1));
      poisson.setLevelBC(lev, phi[lev].get());
      laps.emplace_back(grids[lev], dmap[lev], 1, 1);
    }
    MLMG mlmg(poisson);
    mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
    mlmg.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);
    phi.clear();
    for (int lev = 0; lev < Nlev; ++lev) {
      MultiFab gradAlias(gradMF[lev], amrex::make_alias, 0, AMREX_SPACEDIM); //put the gradient in here
      average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
      gradAlias.mult(-1.0);
    }
    //}
    Print() << "Calculating grad"+inNames[idTlocal] << std::endl;
    for (int lev = 0; lev < Nlev; ++lev) {
      phi.push_back(std::make_unique<MultiFab>(indata[lev],amrex::make_alias,idTlocal,1));
      poisson.setLevelBC(lev, phi[lev].get());
    }
    MLMG mlmg2(poisson); 
    mlmg2.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
    mlmg2.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);
    phi.clear();
    for (int lev = 0; lev < Nlev; ++lev) {
      MultiFab gradAlias(gradMF[lev], amrex::make_alias, AMREX_SPACEDIM, AMREX_SPACEDIM); //put the gradient in here
      average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
      gradAlias.mult(-1.0);
    }
      
    Print() << "Calculating W gradient..." << std::endl;
    for (int lev = 0; lev < Nlev; ++lev) {
      phi.push_back(std::make_unique<MultiFab>(deriveMF[lev],amrex::make_alias,0,1));
      poisson.setLevelBC(lev, phi[lev].get());
    }
    MLMG mlmg3(poisson);
    mlmg3.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
    mlmg3.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);
    phi.clear();
    for (int lev = 0; lev < Nlev; ++lev) {
      MultiFab gradAlias(gradMF[lev], amrex::make_alias, 2*AMREX_SPACEDIM, AMREX_SPACEDIM);
      average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
      gradAlias.mult(-1.0);
    }
  
    //so we have the coefficients and gradients in all directions, lets combine into the fluxes for all three directions.
    //Store in derive MF
    for (int lev=0; lev<Nlev; ++lev)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (amrex::MFIter mfi(indata[lev], amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {

        const Box& bx = mfi.tilebox();
	Array4<Real> const& spec_coeff_a = coeffMF[lev].array(mfi,0);
	Array4<Real> const& soret_coeff_a = coeffMF[lev].array(mfi,1);
	Array4<Real> const& molar_coeff_a = coeffMF[lev].array(mfi,2);
	Array4<Real> const& gradY = gradMF[lev].array(mfi,0);
	Array4<Real> const& gradT = gradMF[lev].array(mfi,AMREX_SPACEDIM);
	Array4<Real> const& gradW = gradMF[lev].array(mfi,2*AMREX_SPACEDIM);

	Array4<Real> const& flux_a = deriveMF[lev].array(mfi,1);
	
	amrex::ParallelFor(bx, [=]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
	  for (int d = 0; d < AMREX_SPACEDIM; d++) {
	    flux_a(i,j,k,d) = spec_coeff_a(i,j,k)*gradY(i,j,k,d);
	    flux_a(i,j,k,d+AMREX_SPACEDIM) = soret_coeff_a(i,j,k)*gradT(i,j,k,d);
	    flux_a(i,j,k,d+2*AMREX_SPACEDIM) = molar_coeff_a(i,j,k)*gradW(i,j,k,d);
	  }
	});

			   
    
      }
    }

    //we now have the diffusive fluxes, gradients and the chemical source term. Let's compute divergence of diffusive flux...
    for (int n =0; n<AMREX_SPACEDIM; n++) {
      Print() << "Calculating gradFspecies"+std::to_string(n) << std::endl;
      for (int lev = 0; lev < Nlev; ++lev) {
	phi.push_back(std::make_unique<MultiFab>(deriveMF[lev],amrex::make_alias,n+1,1)); //plus 1 because 1 is W
	poisson.setLevelBC(lev, phi[lev].get());
	laps.emplace_back(grids[lev], dmap[lev], 1, 1);
      }
      MLMG mlmg(poisson);
      mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
      mlmg.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);
      phi.clear();
      for (int lev = 0; lev < Nlev; ++lev) {
	MultiFab gradAlias(grad2MF[lev], amrex::make_alias, n*AMREX_SPACEDIM, AMREX_SPACEDIM); //put the gradient in here
	average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
	gradAlias.mult(-1.0);
      }
    }
    for (int n =0; n<AMREX_SPACEDIM; n++) {
      Print() << "Calculating gradFsoret"+std::to_string(n) << std::endl;
      for (int lev = 0; lev < Nlev; ++lev) {
	phi.push_back(std::make_unique<MultiFab>(deriveMF[lev],amrex::make_alias,n+1+AMREX_SPACEDIM,1)); //plus 1 because 1 is W
	poisson.setLevelBC(lev, phi[lev].get());
	laps.emplace_back(grids[lev], dmap[lev], 1, 1);
      }
      MLMG mlmg(poisson);
      mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
      mlmg.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);
      phi.clear();
      for (int lev = 0; lev < Nlev; ++lev) {
	MultiFab gradAlias(grad2MF[lev], amrex::make_alias, n*AMREX_SPACEDIM + AMREX_SPACEDIM*AMREX_SPACEDIM, AMREX_SPACEDIM); //put the gradient in here
	average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
	gradAlias.mult(-1.0);
      }
    }
    for (int n =0; n<AMREX_SPACEDIM; n++) {
      Print() << "Calculating gradFmolar"+std::to_string(n) << std::endl;
      for (int lev = 0; lev < Nlev; ++lev) {
	phi.push_back(std::make_unique<MultiFab>(deriveMF[lev],amrex::make_alias,n+1+2*AMREX_SPACEDIM,1)); //plus 1 because 1 is W
	poisson.setLevelBC(lev, phi[lev].get());
	laps.emplace_back(grids[lev], dmap[lev], 1, 1);
      }
      MLMG mlmg(poisson);
      mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
      mlmg.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);
      phi.clear();
      for (int lev = 0; lev < Nlev; ++lev) {
	MultiFab gradAlias(grad2MF[lev], amrex::make_alias, n*AMREX_SPACEDIM + 2*AMREX_SPACEDIM*AMREX_SPACEDIM, AMREX_SPACEDIM); //put the gradient in here
	average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
	gradAlias.mult(-1.0);
      }
    }


    //now we have all the second derivatives, lets combine everything for Sd
    for (int lev=0; lev<Nlev; ++lev)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (amrex::MFIter mfi(indata[lev], amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {

        const Box& bx = mfi.tilebox();
	Array4<Real> const& sd_a = outdata[lev].array(mfi,0);

	Array4<Real> const& gradY_a = gradMF[lev].array(mfi,0);
	Array4<Real> const& diff_spec_a = grad2MF[lev].array(mfi,0);
	Array4<Real> const& diff_soret_a = grad2MF[lev].array(mfi,AMREX_SPACEDIM*AMREX_SPACEDIM);
	Array4<Real> const& diff_molar_a = grad2MF[lev].array(mfi,2*AMREX_SPACEDIM*AMREX_SPACEDIM);
	
	Array4<Real> const& c_a = indata[lev].array(mfi,idClocal);
	Array4<Real> const& rho_a = indata[lev].array(mfi,idRlocal);

	amrex::ParallelFor(bx, [=]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
	  Real maggrad = std::sqrt(AMREX_D_TERM(gradY_a(i,j,k,0)*gradY_a(i,j,k,0),+gradY_a(i,j,k,1)*gradY_a(i,j,k,1),+gradY_a(i,j,k,2)*gradY_a(i,j,k,2)));	  
	  Real lapdiff_spec = AMREX_D_TERM(diff_spec_a(i,j,k,0), + diff_spec_a(i,j,k,AMREX_SPACEDIM + 1), + diff_spec_a(i,j,k,2*AMREX_SPACEDIM+2));
	  Real lapdiff_soret = AMREX_D_TERM(diff_soret_a(i,j,k,0), + diff_soret_a(i,j,k,AMREX_SPACEDIM + 1), + diff_soret_a(i,j,k,2*AMREX_SPACEDIM+2));
	  Real lapdiff_molar = AMREX_D_TERM(diff_molar_a(i,j,k,0), + diff_molar_a(i,j,k,AMREX_SPACEDIM + 1), + diff_molar_a(i,j,k,2*AMREX_SPACEDIM+2));
	  
	  if (maggrad > 1e-4) {
	    sd_a(i,j,k,0) = -lapdiff_spec/(rho_a(i,j,k)*maggrad);
	    sd_a(i,j,k,1) = -lapdiff_soret/(rho_a(i,j,k)*maggrad); 
	    sd_a(i,j,k,2) = -lapdiff_molar/(rho_a(i,j,k)*maggrad);
	    sd_a(i,j,k,3) = -c_a(i,j,k)/(rho_a(i,j,k)*maggrad);
	    sd_a(i,j,k,4) = sd_a(i,j,k,0) + sd_a(i,j,k,1) + sd_a(i,j,k,2) + sd_a(i,j,k,3);
	    sd_a(i,j,k,5) = sd_a(i,j,k,4)*rho_a(i,j,k)/rhou;
	  } else {
	    sd_a(i,j,k,0) = 0.0;
	    sd_a(i,j,k,1) = 0.0;
	    sd_a(i,j,k,2) = 0.0;
	    sd_a(i,j,k,3) = 0.0;
	    sd_a(i,j,k,4) = 0.0;
	    sd_a(i,j,k,5) = 0.0;
	  }
	});


      }
    }
    

    
    std::string outfile(getFileRoot(plotFileName) + "_Sd");
    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(outfile, Nlev, GetVecOfConstPtrs(outdata), outNames,
                                   geoms, 0.0, isteps, refRatios);
  }
  Finalize();
  return 0;
}
