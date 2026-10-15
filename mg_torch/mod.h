#pragma once

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

struct phi_init
{
  double phi_init = 0.0;
  bool got = false;
  bool vanish(double phi) {
    if (phi_init > 0.0)
      return phi <= 0.0;
    else
      return phi >= 0.0;
  }
};

/* Data at end points */
struct EndPoint
{
  double p_son[MAXQ]; /* Sonic end point */
  double p_phi[MAXQ]; /* Phi=0 end point */
  double dp; /*det speed*/
  /* Boundary det speeds */
  double m_dp;
  double phi_dp;
  /* Mach numbers */
  double m_son;
  double m_phi;
  /* Phi */
  double phi_son;
  double phi_phi;
  double xloc;
};

enum eSol
{
  Bisect = 0,
  Single = 1
};

struct SolParams
{
  double k = 0.0; // No curvature
  bool helium = true;
  double dn = 1.4e9;
  double dx = 1e-16; // Initial step size
  double nsteps = -1;
  double target = -1; // End target
  std::string react = "ISO13";
  std::string exp = "timmesHydrostatic";
  double rho0 = 1e9;
  double T0 = 5e9;
  std::string init = "TIMMES";

  void fromFile(const std::string& fileName)
  {
    std::ifstream ifs(fileName.c_str());
    std::string line;
    std::vector<std::string> tokens;

    while (std::getline(ifs, line)) {
      std::stringstream ss; ss << line;
      std::string out;
      tokens.clear();
      while (std::getline(ss,out, ' ')) {
        tokens.emplace_back(out);
      }
      auto cmd = tokens[0];
      if (cmd == "react") {
        react = tokens[1];
      } 
      else if (cmd == "exp") {
        exp = tokens[1];
      } 
      else if (cmd == "dn") {
        dn = atof(tokens[1].c_str());
      } 
      else if (cmd == "k") {
        k = atof(tokens[1].c_str());
      } 
      else if (cmd == "dx") {
        dx = atof(tokens[1].c_str());
      } 
      else if (cmd == "target") {
        target = atof(tokens[1].c_str());
      } 
      else if (cmd == "rho0") {
        rho0 = atof(tokens[1].c_str());
      } 
      else if (cmd == "T0") {
        T0 = atof(tokens[1].c_str());
      }
      else if (cmd == "init") {
        init = tokens[1];
      }
      else {
        std::cout << "Unknown cmd " << cmd << std::endl;
        abort();
      }
    }
  }

  void fromArgs(int argc,char* argv[]) {
    if (argc > 2) {
      k = atof(argv[2]);
    }
    if (argc > 3) {
      helium = atoi(argv[3]);
    }
    if (argc > 4) {
      dn = atof(argv[4]);
    }
    if (argc > 5) {
      dx = atof(argv[5]);
    }
    if (argc > 6) {
      nsteps = atoi(argv[6]);
    }
  }

  void Out()
  {
    std::cout << " --- ARGS ---- " << std::endl;
    std::cout << "k: " << k << " he? " << helium << " dn: " << dn << " dx: " << dx << " nsteps: " << nsteps << std::endl;
    std::cout << " ------------- " << std::endl;
  }
};

class Timer
{
public:
  using T_Time = std::chrono::time_point<std::chrono::system_clock>;
  T_Time m_start;
  T_Time m_end;

  double elapsed()
  {
    return std::chrono::duration_cast<std::chrono::microseconds>(m_end - m_start).count();
  }
  void start()
  {
    m_start = std::chrono::system_clock::now();
  }
  void stop()
  {
    m_end = std::chrono::system_clock::now();
  }

};
