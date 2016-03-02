/*
 * multiplane.cpp
 *
 *  Created on: Jan 11, 2012
 *      Author: bmetcalf
 *
 *      This program is for testing the adaptive griding on a random field by adapting to the
 *      high convergence regions.
 */

#include <slsimlib.h>
#include <sstream>
#include <string.h>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <omp.h>
#include <thread>
#include <mutex>

using namespace std;

static std::mutex barrier;

int main(int arg,char **argv){
  
  bool do_caustics = true,do_wholeimage = false,do_stamps = false,do_qso = false;
  
  //const double gal_range = 1.10*pi*sqrt(source.getFOV())/180; //2*0.18*pi/180;
  const double gal_range = 0.00383916446299277;
  const double resolution = 0.1*pi/180/60/60;  // Euclid VIS
  //const double resolution = 0.2*pi/180/60/60;  // KIDS pixels
  //const double resolution = 0.3*pi/180/60/60;    // Euclid IR pixels
  const int Nsourcesperlens= 3;                // number os sources implanted per lens
  const int Nsubfields = 8;                       // the total region is broken into Nfields*Nfields subfields for caustic finding
  //const double source_dist_radius = 0.75;   // for the stamps the sources are distributed within this number times the caustis size
  const double source_dist_radius = 1.5;   // for the stamps the sources are distributed within this number times the caustis size
  const double sourceMag_limit_Iband = 26; //
  //long noise_seed = -1*time(NULL);
  long noise_seed = -1234567+22; // +22;
  
  
  //unsigned int
  long seed = 12213-22;
  Utilities::RandomNumbers_NR random(seed);
  
  
  const bool make_noise = true;
  Telescope telescope = Euclid_VIS;
  //Telescope telescope = Euclid_J;
  //Telescope telescope = Euclid_H;
  
  Observation observe(telescope);
  
  //   make initial grid
  double center[2] = {0.0,0.0};
  char chrKcaustic[100];
  
  printf("initializing model\n");
  string paramfile;
  
  if(arg > 1) paramfile.assign(argv[1],strlen(argv[1]));
  else paramfile = "ParamFiles/param";
  cout << "using parameter file: " << paramfile << endl;
  
  double z_caustics = 3.0;         // redshift of sources for caustic finding
  if(arg > 2) sscanf(argv[2],"%lf",&z_caustics);
  
  cout << "Create model" << endl;
  
  InputParams params(paramfile);
  //SourceMultiAnaGalaxy source(params);
  
  Lens lens(params,&noise_seed);
  // Set catalog files
  string output,band_str;
  Band band;
  
  // ************* create log file ******************************
  
  params.get("outputfile",output);
  ofstream logfile;
  time_t now = time(0);
  struct tm date = *localtime(&now);
  snprintf(chrKcaustic,100,"%f",z_caustics);
  logfile.open (output + "_zcaust" + chrKcaustic + ".log");
  logfile << "Log file, multi_plane_field.cpp \n";
  logfile << date.tm_hour << ":" << date.tm_min << "   " << date.tm_mday << "/" << date.tm_mon << "/" << date.tm_year << endl;
  logfile << "do_caustics: " << do_caustics << endl;
  logfile << "do_wholeimage: " << do_wholeimage << endl;
  logfile << "do_stamps: " << do_stamps << endl;
  logfile << "do_qso: " << do_qso << endl;
  logfile << "parameter file: " << params.filename() << endl;
  logfile << "gal_range: " << gal_range << " radians" << endl;
  logfile << "resolution: " << resolution*180*60*60/pi << " arcsec" << endl;
  logfile << "z_caustics: " << z_caustics << endl;
  logfile << "Nsourcesperlens: " << Nsourcesperlens << endl;
  logfile << "Nfields :" <<  Nsubfields << endl;
  logfile << "source_dist_radius: " << source_dist_radius << endl;
  logfile << "sourceMag_limit_Iband: " << sourceMag_limit_Iband << endl;
  logfile << "make_noise: " << make_noise << endl;
  if(make_noise){
    logfile << "telescope " << telescope << endl;
    logfile << "noise_seed: " << noise_seed << endl;
  }
  logfile << "seed for RandomNumbers: " << seed << endl;
  
  //*****************************************
  
  
  string tmp_string ;
  
  params.get("source_band",band_str);
  params.get("source_band",band);
  
  output = output + band_str;
  
  time_t t0,t1,t2;
  time(&t0);
  
  params.print_used();
  // cout << source.getZ() << endl;
  
  
  //cout << "Field of view from source = " << source.getFOV() << " from gal_rnage " << pow(gal_range*180/pi,2) << " deg^2 "<< endl;
  //  logfile << "Field of view from source = " << source.getFOV() << " from gal_range " << pow(gal_range*180/pi,2) << " deg^2 "<< endl;
  
  //const double sim_range = pi*sqrt(lens->getfov()/pi)*2/180; // convert to radians
  //cout << "Range of lensing simulation box is " << 180*sim_range/pi << " degrees, area is "
  //		  << lens->getfov() << " square degrees" << endl;
  //double sim_range = pi*sqrt(lens.getfov()/pi)*2/180;// // convert to radians
  double boxsize= 9.3e-06;
  //double sim_range = pi* boxsize*2/180;// // convert to radians
  
  double sim_range = boxsize*2.*10;
  
  //if(sim_range == 0.0) sim_range = sqrt(source.getFOV()/pi)*pi/180;
  
  cout << "Range of lensing simulation box is " << 180*sim_range/pi << " degrees, area is "
  << lens.getfov() << " square degrees" << endl;
  logfile << "Range of lensing simulation box is " << 180*sim_range/pi << " degrees, area is "
  << lens.getfov() << " square degrees" << endl;
  
  //  params->print_unused();
  
  //range = 3.0*pi/180/60; // range of image
  //long Npixels = prevpower((long)(range/res + 1));
  //long Npixels = (long)pow((float)2,12);
  int Ncaustics,j;
  bool dummy;
  
  
  //range = 0.18*resolution*Npixels;
  
  cout << "range of view for galaxies is " << gal_range*180/pi << " deg" << endl;
  logfile << "range of view for galaxies is " << gal_range*180/pi << " deg" << endl;
  // initialize tree and lens model
  
  long Npixels;
  //CausticData causticdata(0);
  
  // TODO: test line
  //***  find all caustics
  
  
  std::vector<ImageFinding::CriticalCurve> caustics;
  CausticDataStore cdatastore(caustics);
  
  int Ncrit = 0,Ncrit_tmp;
  PosType rmax,rmin,rave;
  
  
  lens.ResetSourcePlane(z_caustics,true);
  
  Npixels = 256; //38; // Utilities::prevpower((long)(sim_range/resolution + 1))/2/2;///2;// /2;///2;
  // Npixels = Utilities::prevpower((long)(sim_range/resolution + 1))/2/2/2 /2;///2;
  logfile << "finding all critical curves above scale " << 180*60*60*sim_range/Npixels/pi << " arcsec"<< endl;
  cout << "    initializing Grid" << endl;
  
  //catalog_caustic << "# " << " all critical lines above a scale of " << 180*60*60*sim_range/Npixels/pi << " arcsec,  field of view: " << lens->fieldofview << " square degrees" << std::endl;
  
  time(&t1);
  cout << "   finding caustics" << endl;
  
  double tmp_center[2],tmp_center_map[2],maxEradius=0;
  // int jold=0;
  Ncaustics = Ncrit = j = 0;
  
  // int i_max=0;
  
  double pix_size = 1.*.1/60./60./180.*pi;
  
  //tmp_center[0] = -0.00163436; // *pi/180; //-0.00161226*pi/180.; // from deg in rad
  //tmp_center[1] = -0.00012315; //*pi/180.; // -0.000124854*pi/180.;
  //tmp_center[0] =  -0.00153108; // -0.0016125// <- coord of biggest in smallfield ;
  //tmp_center[1] = -0.0014296; // -0.000125061  //<- coord of biggest in smallfield ;
  
  tmp_center[0] = -0.00015583; // largest area
  tmp_center[1] = 0.00167417;
  tmp_center_map[0] =  -0.000157517 ; // caust center coords
  tmp_center_map[1] = 0.00164907-0.00005;
  
  //tmp_center[0] =-3.60571e-05;
  //tmp_center[1] =0.000477221 ; // radial orphan
  PixelMap map(tmp_center,Npixels,sim_range/Npixels);
  PixelMap mapc(tmp_center_map,Npixels,sim_range/2./Npixels);
  //Utilities::PositionFromIndex(i_subfield,tmp_center,Nsubfields,sim_range*(1 - 1.0/Nsubfields) ,center);
  
  std::cout << " center = " << tmp_center[0] << "  " << tmp_center[1] << endl;
  std::cout << " range is " << sim_range << " resolution (asec) = " << 180*60*60*sim_range/Npixels/pi  << endl;
  
  //if(i_subfield > 0) delete grid;
  Grid grid(&lens,Npixels,tmp_center,sim_range);
  
  
  
  //******** test lines
  //grid.ReShoot(&lens);
  //grid.writeFits(tmp_center, Npixels/Nfields, sim_range*Nfields/Npixels, kappa, "test_rand");
  //grid.writeFits(tmp_center, Npixels/Nfields, sim_range*Nfields/Npixels, invmag, "test_rand");
  //exit(1);
  // *********************************************/
  
  std::vector<ImageFinding::CriticalCurve> tmp_caustics;
  
  
  //ImageFinding::find_crit(&lens,&grid,tmp_caustics,&Ncrit_tmp,0.03*pi/60/60/180,&dummy,true,true,1.0/mu_min);
  ImageFinding::find_crit(&lens,&grid,tmp_caustics,&Ncrit_tmp,0.03*pi/60/60/180,0,false);
  
  
  
  // find largerst Einstein radius in this batch
  for(size_t ii=0;ii<tmp_caustics.size();++ii){
    //std::cout<<   "HERE "<< std::endl;
    tmp_caustics[ii].z_source = z_caustics;
    tmp_caustics[ii].CriticalRadius(rmax,rmin,rave);
    mapc.AddCurve(tmp_caustics[ii].caustic_curve_outline,1);
    map.AddCurve(tmp_caustics[ii].critical_curve,1);
    maxEradius = MAX(maxEradius,rave);
  }
  
  
  
  
  /// copy caustic data to output storage structure
  cdatastore.addcrits(tmp_caustics);
  
  /// add new caustics to list of old ones
  caustics.insert(caustics.end(), tmp_caustics.begin(),tmp_caustics.end());
  
  Ncaustics = Ncrit = caustics.size();
  //find_crit2(&lens,&grid,caustics,&Ncrit_tmp,0.03*pi/60/60/180,&dummy,true,true,0.0);
  //Ncrit += Ncrit_tmp;
  
  time(&t2);
  cout << endl << "found "  << Ncrit << " critical curves and " << Ncaustics << " caustic curves in " << difftime(t2,t1) << " sec" << endl;
  logfile << endl << "found "  << Ncrit << " critical curves and " << Ncaustics << " caustic curves in " << difftime(t2,t1) << " sec" << endl;
  
  //*********************/
  //jold = causticdata.numberOfCaustics();
  //causticdata.resize(Ncaustics);
  
  tmp_string = output + "_caustic_catalog_z" + std::to_string(z_caustics)
  + "_sub";
  cout << "  printing " << tmp_string << endl;
  
  // log results
  cdatastore.printfile(tmp_string,paramfile,lens.getfov(),sim_range/Npixels);
  std::cout << "FOV: " << lens.getfov() <<  "sim_range/Npixels: " << sim_range/Npixels  <<  std::endl;
  map.printFITS("!max_crit.fits");
  mapc.printFITS("!max_caust.fits");
  grid.writeFits(tmp_center,Npixels,pix_size,KAPPA,"!test");
  grid.writeFits(tmp_center,Npixels,pix_size,INVMAG,"!test");
  grid.writeFits(tmp_center,Npixels,pix_size,ALPHA1,"!test");
  
  
  
  params.get("outputfile",output);
  tmp_string = output + "_caustic_catalog_z" + std::to_string(z_caustics);
  cout << "  printing " << tmp_string << endl;
  cdatastore.printfile(tmp_string,paramfile,lens.getfov(),sim_range/Npixels);
  
  
  cout << "Largest Einstein radius " << maxEradius*180*60*60/pi << "  arcsec" << endl;
  //  cout << "Number of galaxies: " << source.getNumberOfGalaxies() << endl;
  logfile << "Largest Einstein radius " << maxEradius*180*60*60/pi << "  arcsec" << endl;
  // logfile << "Number of galaxies: " << source.getNumberOfGalaxies() << endl;
  
  //cout << lens.getmaxalpha() << endl;
 // cout << "Done in " << std::setprecision(2) << difftime(t1,t0)/60. << " mins" << "\n\n" << endl;
  

  return 0;
}


