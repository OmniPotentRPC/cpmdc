#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
extern "C" double native_step(double,double,double);
int main() {
  try {
    for(int exponent=-300;exponent<=300;exponent+=5) {
      const double scale=std::pow(10.,exponent);
      for(double direction:{-1.,1.})for(double trial:{.125,1.,3.}) {
        const double actual=native_step(2*direction*scale,-direction*scale,trial);
        const double expected=2*trial/3;
        if(!std::isfinite(actual)||std::abs(actual-expected)>2e-15*std::max(1.,expected))
          throw std::runtime_error("Native line-search root depends on derivative scale");
      }
    }
    if(native_step(0,0,1)!=1||native_step(1,1,.5)!=.5)
      throw std::runtime_error("Degenerate derivative pair must retain the trial step");
    if(std::abs(native_step(2,1,1)-2)>2e-15||std::abs(native_step(1,2,1)+1)>2e-15)
      throw std::runtime_error("Native extrapolation sign or magnitude differs");
    std::cout<<"726 scaled native line-search roots and four degenerate/extrapolation cases pass\n";
    return 0;
  }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
