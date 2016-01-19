#include "CMGTools/TTHAnalysis/interface/DistributionRemapper.h"

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <TH1.h>
#include <TAxis.h>

DistributionRemapper::DistributionRemapper(const TH1 *source, const TH1 *target, bool bounded, bool twoWay, bool spline, bool debug) :
    xmin_(source->GetXaxis()->GetXmin()),
    ymin_(target->GetXaxis()->GetXmin()),
    xmax_(source->GetXaxis()->GetXmax()),
    ymax_(target->GetXaxis()->GetXmax()),
    x_(source->GetNbinsX()+1),
    y_(source->GetNbinsX()+1),
    linear_(!spline),
    interp_(0)
{
    int ns = source->GetNbinsX();
    int nt = target->GetNbinsX(); 
    const TAxis *axs = source->GetXaxis();
    const TAxis *axt = target->GetXaxis(); 

    std::vector<double> xt, yt;
    makeCdf(xt, yt, target, bounded);
    if (debug) {
        printf("target cdf: \n");
        for (unsigned int i = 0; i < xt.size(); ++i) {
            printf("xt yt %3d % 10.5f % 10.5f\n", i, xt[i], yt[i]);
        }
    }
    ROOT::Math::Interpolator tinv(yt, xt, ROOT::Math::Interpolation::kLINEAR);
   
    double tnorm = 1.0/target->Integral(0, nt+1);
    double snorm = 1.0/source->Integral(0, ns+1);

    // loop over bin
    //printf("preparing dataset for morphing function (%d bins)\n", ns);
    double srun = bounded ? 0 : snorm*source->GetBinContent(0);
    x_[0] = axs->GetBinLowEdge(1);
    y_[0] = bounded ? 0 : tinv.Eval(srun);
    for (int i = 1; i <= ns; ++i) {
        srun += tnorm * source->GetBinContent(i);
        x_[i] = axt->GetBinUpEdge(i);
        y_[i] = tinv.Eval(std::min(1.,srun));
    }
    if (twoWay) { 
        if (debug) {
            printf("one-way mapping: \n");
            for (unsigned int i = 0; i < x_.size(); ++i) {
               printf("x  y  %3d % 10.5f % 10.5f\n", i, x_[i], y_[i]);
            }
        }
        std::vector<std::pair<double,double> > v(x_.size());
        for (unsigned int i = 0; i < x_.size(); ++i) {
            v[i] = std::make_pair(x_[i], y_[i]);
        }
        // now add additional points for the CDF of the source
        std::vector<double> xs, ys;
        makeCdf(xs, ys, source, bounded);
        if (debug) {
            printf("source cdf: \n");
            for (unsigned int i = 0; i < xs.size(); ++i) {
                printf("xs ys %3d % 10.5f % 10.5f\n", i, xs[i], ys[i]);
            }
        }
        ROOT::Math::Interpolator sinv(ys, xs, ROOT::Math::Interpolation::kLINEAR);
        double trun = bounded ? 0 : tnorm*target->GetBinContent(0);
        double x = axs->GetBinLowEdge(1);
        double y = bounded ? 0 : sinv.Eval(srun);
        if (std::find(x_.begin(), x_.end(), y) == x_.end()) {
            v.emplace_back(y,x);
        }
        for (int i = 1; i < ns; ++i) {
            trun += snorm * target->GetBinContent(i);
            x = axs->GetBinUpEdge(i);
            y = sinv.Eval(std::min(1.,trun));
            if (std::find(x_.begin(), x_.end(), y) == x_.end()) {
                v.emplace_back(y,x);
            }
        }
        std::sort(v.begin(), v.end());
        x_.resize(v.size());
        y_.resize(v.size());
        for (unsigned int i = 0; i < x_.size(); ++i) {
            x_[i] = v[i].first;
            y_[i] = v[i].second;
        }
    }
    if (debug) {
        printf("final mapping: \n");
        for (unsigned int i = 0; i < x_.size(); ++i) {
           printf("x  y  %3d % 10.5f % 10.5f\n", i, x_[i], y_[i]);
        }
    }
}

void DistributionRemapper::makeCdf(std::vector<double> &xt, std::vector<double> &yt, const TH1* target, bool bounded) 
{
    int nt = target->GetNbinsX(); 
    const TAxis *axt = target->GetXaxis(); 
    xt.resize(nt+3);
    yt.resize(nt+3);

    // make inverse cdf of target
    //printf("preparing dataset for inverse cdf of target (%d bins)\n", nt);
    double tnorm = 1.0/target->Integral(0, nt+1);
    xt[0] = axt->GetXmin() - (!bounded ? 0.5*(axt->GetXmax()-axt->GetXmin()) : 0);
    yt[0] = 0.0;
    int j = 0;
    if (!bounded && target->GetBinContent(0) != 0) {
        j = 1;
        xt[1] = axt->GetBinLowEdge(1);
        yt[1] = tnorm*target->GetBinContent(0);
    }
    for (int i = 1; i <= (bounded ? nt - 1 : nt); ++i) {
        if (target->GetBinContent(i) == 0) continue;
        j++;
        xt[j] = axt->GetBinUpEdge(i);
        yt[j] = yt[j-1] + tnorm*target->GetBinContent(i);
        if (yt[j] >= 1) break;
    }
    if (yt[j] < 1) {
        j++;
        yt[j] = 1.0;
        xt[j] = axt->GetXmax() + (!bounded ? 0.5*(axt->GetXmax()-axt->GetXmin()) : 0);
    }
    xt.resize(j+1);
    yt.resize(j+1);
}


DistributionRemapper::~DistributionRemapper() 
{
    delete interp_; interp_ = 0;
}

double DistributionRemapper::Eval(double x) const 
{
    if (!interp_) init();
    if (x < xmin_) return ymin_;
    if (x > xmax_) return ymax_;
    return interp_->Eval(x);
}

void DistributionRemapper::init() const 
{
    interp_ = new ROOT::Math::Interpolator(x_, y_, linear_ ? ROOT::Math::Interpolation::kLINEAR : ROOT::Math::Interpolation::kCSPLINE);
}


