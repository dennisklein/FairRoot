// ProgressLog.h is a part of the PYTHIA event generator.
// Copyright (C) 2026 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#ifndef PROGRESSLOG_H
#define PROGRESSLOG_H

#include "Pythia8/Basics.h"
#include <unistd.h>
#include <ctime>
#include <string>
#include <iostream>
#include <iomanip>

namespace Pythia8 {

using namespace std;

// ProgressLog is a simple helper class to monitor the progress of a run.
// When used in the main event loop, it will with suitably (logarithmically)
// spaced intervals, print out one line with information about the number
// of events generated, two estimates (based on instantaneous and average
// CPU consumption) of when the run will be completed, the host on which
// the program is run and its process number.

class ProgressLog {

public:

  // Create an object for monitoring the progress of a run with NIn iterations.
  ProgressLog( long NIn, int maxsec = 1200) : secstep(maxsec) {
    init(NIn);
  }

  // Intermittently write out a line of progress information, giving
  // the current iteration (in the interval [0:N[ ).
  void operator()(int cnt) {
    tick(cnt + 1, N);
  }

  // Intermittently write out a line of progress information using the
  // internal counter.
  void operator()() {
    tick(++count, N);
  }

  // Intermittently write out a line of progress information giving
  // the current iteration and the total number of iterations.
  void tick(long i, long n) {
    if ( !statusTime(i, n) ) return;

    // Wall-clock and CPU time in seconds since the start, and the date now.
    double walli = wallTimer.elapsed() / 1000.;
    double cpui = cpuTimer.elapsed() / 1000.;
    time_t timei = time(0);

    // CPU efficiency since the last status line and since the start.
    double eff = 1.0;
    if ( walli > wall1 && cpui > cpu1 )
      eff = (cpui-cpu1)/(walli-wall1);
    if ( eff >= 1.0 ) eff = 0.999999;
    int ieff = 100*eff;
    double eff0 = 1.0;
    if ( walli > 0. && cpui > 0. )
      eff0 = cpui/walli;
    if ( eff0 >= 1.0 ) eff0 = 0.999999;
    int ieff0 = 100*eff0;

    // Estimated remaining CPU time and expected time of completion.
    double cpuRem = (cpui*(n-i))/i;
    time_t timen = (time_t)(timei + cpuRem/eff + 30.0);
    time_t timen0 = (time_t)(timei + cpuRem/eff0 + 30.0);
    char date[1024];
    char daten[1024];
    char daten0[1024];
    strftime(date,1024,"%y.%m.%d %H:%M",localtime(&timei));
    strftime(daten,1024,"%H:%M",localtime(&timen));
    strftime(daten0,1024,"%H:%M",localtime(&timen0));
    long ii = i;
    if ( n - i < n/10 ) ii = i - n;
    time_t dayn = (timen - timei)/86400;
    time_t dayn0 = (timen0 - timei)/86400;

    ostream & os = cout;

    if ( dayn <= 0 && dayn0 <= 0 ) {
      os << date << " " << setw(8) << ii << "/" << setw(9);
      os.setf(ios::left, ios::adjustfield);
      os << n << " etc:   " << daten << "[";
      os.setf(ios::right, ios::adjustfield);
      os << setw(2) << ieff << "%]   " << daten0 << "[" << ieff0 << "%] "
         << host << ":" << pid << endl << flush;
    } else {
      os << date << " " << setw(8) << ii << "/" << setw(9);
      os.setf(ios::left, ios::adjustfield);
      os << n << " etc: " << dayn << "+" << daten << "[";
      os.setf(ios::right, ios::adjustfield);
      os << setw(2) << ieff << "%] "
         << dayn0 << "+" << daten0 << "[" << ieff0 << "%] "
         << host << ":" << pid << endl << flush;
    }

    cpu1 = cpui;
    wall1 = walli;

  }

  // Check if this is a good time to print out a status line.
  bool statusTime(long i, long n) const {
    if ( i <= 0 ) return false;
    if ( i == n ) return true;
    if ( i > n/2 ) i = n-i;
    while ( i >= 10 && !(i%10) ) i /= 10;
    if ( i == 1 || i == 2 || i == 5 ) return true;
    return secstep > 0 && wallTimer.elapsed()/1000. > wall1 + secstep;
  }

  // Initialise the basic engine.
  void init(long n) {
    N = n;
    count = 0;
    wallTimer.start();
    cpuTimer.start();
    wall1 = cpu1 = 0.;
    time_t time0 = time(0);
    char name[1024];
    gethostname(name,1024);
    host = name;
    if ( host.find(".") != string::npos )
      host = host.substr(0, host.find("."));
    pid = getpid();
    char date[1024];
    strftime(date,1024,"%y.%m.%d %H:%M",localtime(&time0));
    ostream & os = cout;
    os << date << "        0/" << setw(9);
    os.setf(ios::left, ios::adjustfield);
    os << n;
    os.setf(ios::right, ios::adjustfield);
    os << " Initializing...                "
       << host << ":" << pid << endl << flush;
  }

private:

  // If larger than 0, a status line will be written every secstep
  // second.
  int secstep;

  // Timers for the wall-clock and CPU time since the run was started.
  Timer wallTimer{Timer::WALL}, cpuTimer{Timer::CPU};

  // The wall-clock and CPU time, in seconds since last status line.
  double wall1, cpu1;

  // The host on which we are running.
  string host;

  // The pid of the current process.
  pid_t pid;

  // The number of iterations
  long N;

  // The number of iterations so far
  long count;

};

}

#endif
