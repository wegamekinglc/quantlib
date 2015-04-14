#include <ql/quantlib.hpp>

#include <boost/assign/std/vector.hpp>
#include <boost/make_shared.hpp>
#include <boost/timer.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

using namespace boost::assign;
using namespace QuantLib;

typedef Real dbl;
typedef CppAD::AD<double> dblAD;

class Timer {
    boost::timer timer_;
    double elapsed_;

  public:
    void start() { timer_ = boost::timer(); }
    void stop() { elapsed_ = timer_.elapsed(); }
    double elapsed() const { return elapsed_; }
};

int main() {

    Timer timer;

    // evaluation date

    Date refDate(13, April, 2015);
    Settings::instance().evaluationDate() = refDate;

    // yts

    timer.start();

    Handle<Quote_t<dbl> > swapQuote(
        boost::make_shared<SimpleQuote_t<dbl> >(0.03));
    Handle<Quote_t<dblAD> > swapQuoteAD(
        boost::make_shared<SimpleQuote_t<dblAD> >(0.03));
    boost::shared_ptr<Euribor_t<dbl> > euribor6mBt =
        boost::make_shared<Euribor_t<dbl> >(6 * Months);
    boost::shared_ptr<Euribor_t<dblAD> > euribor6mBtAD =
        boost::make_shared<Euribor_t<dblAD> >(6 * Months);
    std::vector<boost::shared_ptr<RateHelper_t<dbl>::Type> > helpers;
    std::vector<boost::shared_ptr<RateHelper_t<dblAD>::Type> > helpersAD;
    for (Size i = 1; i <= 30; ++i) {
        boost::shared_ptr<RateHelper_t<dbl>::Type> tmp =
            boost::make_shared<SwapRateHelper_t<dbl> >(
                swapQuote, i * Years, TARGET(), Annual, ModifiedFollowing,
                Thirty360(), euribor6mBt);
        helpers.push_back(tmp);
        boost::shared_ptr<RateHelper_t<dblAD>::Type> tmpAD =
            boost::make_shared<SwapRateHelper_t<dblAD> >(
                swapQuoteAD, i * Years, TARGET(), Annual, ModifiedFollowing,
                Thirty360(), euribor6mBtAD);
        helpersAD.push_back(tmpAD);
    }
    boost::shared_ptr<YieldTermStructure_t<dbl> > yts6m = boost::make_shared<
        PiecewiseYieldCurve<ZeroYield, Linear, IterativeBootstrap, dbl> >(
        refDate, helpers, Actual365Fixed());
    boost::shared_ptr<YieldTermStructure_t<dblAD> > yts6mAD =
        boost::make_shared<
            PiecewiseYieldCurve<ZeroYield, Linear, IterativeBootstrap, dblAD> >(
            refDate, helpersAD, Actual365Fixed());

    Handle<YieldTermStructure_t<dbl> > yts6m_h(yts6m);
    Handle<YieldTermStructure_t<dblAD> > yts6mAD_h(yts6mAD);

    yts6m_h->enableExtrapolation();
    yts6mAD_h->enableExtrapolation();

    std::cout << "test yts = " << yts6m_h->discount(1.0) << std::endl;
    std::cout << "test ytsAD = " << yts6mAD_h->discount(1.0) << std::endl;
    timer.stop();
    std::cout << "timing: " << timer.elapsed() << std::endl;

    // index

    timer.start();

    boost::shared_ptr<Euribor_t<dbl> > euribor6m =
        boost::make_shared<Euribor_t<dbl> >(6 * Months, yts6m_h);
    boost::shared_ptr<Euribor_t<dblAD> > euribor6mAD =
        boost::make_shared<Euribor_t<dblAD> >(6 * Months, yts6mAD_h);

    // SABR swaption volatility structure

    std::vector<Period> optionTenors, swapTenors;

    optionTenors += 3 * Months, 6 * Months, 1 * Years, 2 * Years, 3 * Years,
        4 * Years, 5 * Years, 6 * Years, 7 * Years, 8 * Years, 9 * Years,
        10 * Years, 12 * Years, 15 * Years, 20 * Years, 30 * Years;

    swapTenors += 1 * Years, 2 * Years, 3 * Years, 4 * Years, 5 * Years,
        6 * Years, 7 * Years, 8 * Years, 9 * Years, 10 * Years, 15 * Years,
        20 * Years, 25 * Years, 30 * Years;

    // atm vol structure
    Matrix_t<dbl> atmVols(optionTenors.size(), swapTenors.size(), 0.20);
    boost::shared_ptr<SwaptionVolatilityMatrix_t<dbl> > swatm =
        boost::make_shared<SwaptionVolatilityMatrix_t<dbl> >(
            TARGET(), ModifiedFollowing, optionTenors, swapTenors, atmVols,
            Actual365Fixed());
    Handle<SwaptionVolatilityStructure_t<dbl> > swatm_h(swatm);

    std::cout << "test atm = " << swatm_h->volatility(1.0, 1.0, 0.05)
              << std::endl;
    timer.stop();
    std::cout << "timing: " << timer.elapsed() << std::endl;

    // we do not need spreaded vols, since we assume that
    // we know the sabr parameters

    timer.start();

    std::vector<dbl> strikeSpreads(1, 0.0);
    std::vector<std::vector<Handle<Quote_t<dbl> > > > volSpreads, sabrParams;
    for (Size i = 0; i < optionTenors.size() * swapTenors.size(); ++i) {
        // spreaded vols
        std::vector<Handle<Quote_t<dbl> > > volTmp(
            1, Handle<Quote_t<dbl> >(
                   boost::make_shared<SimpleQuote_t<dbl> >(0.0)));
        volSpreads.push_back(volTmp);
        // sabr parameters
        std::vector<Handle<Quote_t<dbl> > > sabrTmp;
        sabrTmp += Handle<Quote_t<dbl> >(
            boost::make_shared<SimpleQuote_t<dbl> >(0.03)); // alpha
        sabrTmp += Handle<Quote_t<dbl> >(
            boost::make_shared<SimpleQuote_t<dbl> >(0.60)); // beta
        sabrTmp += Handle<Quote_t<dbl> >(
            boost::make_shared<SimpleQuote_t<dbl> >(0.12)); // nu
        sabrTmp += Handle<Quote_t<dbl> >(
            boost::make_shared<SimpleQuote_t<dbl> >(0.30)); // rho
        sabrParams.push_back(sabrTmp);
    }

    std::vector<bool> paramFixed;
    paramFixed += true, true, true, true;

    boost::shared_ptr<SwapIndex_t<dbl> > indexBase =
        boost::make_shared<EuriborSwapIsdaFixA_t<dbl> >(30 * Years, yts6m_h);
    boost::shared_ptr<SwapIndex_t<dbl> > indexBaseShort =
        boost::make_shared<EuriborSwapIsdaFixA_t<dbl> >(2 * Years, yts6m_h);

    boost::shared_ptr<SwaptionVolatilityStructure_t<dbl> > swvol =
        boost::make_shared<SwaptionVolCube1_t<dbl>::Type>(
            swatm_h, optionTenors, swapTenors, strikeSpreads, volSpreads,
            indexBase, indexBaseShort, true, sabrParams, paramFixed, true,
            boost::shared_ptr<EndCriteria_t<dbl> >(), 100.0);
    Handle<SwaptionVolatilityStructure_t<dbl> > swvol_h(swvol);

    std::cout << "test atm = " << swvol_h->volatility(1.0, 1.0, 0.05)
              << std::endl;
    timer.stop();
    std::cout << "timing: " << timer.elapsed() << std::endl;

    // our bermudan to price

    timer.start();

    Date effectiveDate = TARGET().advance(refDate, 2 * Days);
    Date maturityDate = TARGET().advance(effectiveDate, 10 * Years);
    Schedule fixedSchedule(effectiveDate, maturityDate, 1 * Years, TARGET(),
                           ModifiedFollowing, ModifiedFollowing,
                           DateGeneration::Forward, false);
    Schedule floatingSchedule(effectiveDate, maturityDate, 6 * Months, TARGET(),
                              ModifiedFollowing, ModifiedFollowing,
                              DateGeneration::Forward, false);

    dbl strike = 0.03;

    boost::shared_ptr<VanillaSwap_t<dbl> > underlying =
        boost::make_shared<VanillaSwap_t<dbl> >(
            VanillaSwap_t<dbl>::Payer, 1.0, fixedSchedule, strike, Thirty360(),
            floatingSchedule, euribor6m, 0.0, Actual360());
    boost::shared_ptr<VanillaSwap_t<dblAD> > underlyingAD =
        boost::make_shared<VanillaSwap_t<dblAD> >(
            VanillaSwap_t<dblAD>::Payer, 1.0, fixedSchedule, strike,
            Thirty360(), floatingSchedule, euribor6mAD, 0.0, Actual360());

    std::vector<Date> exerciseDates;
    for (Size i = 1; i < 10; ++i) {
        exerciseDates.push_back(TARGET().advance(fixedSchedule[i], -2 * Days));
        // std::cout << "exercise date #" << i << " is " << exerciseDates[i-1]
        // << std::endl;
    }

    boost::shared_ptr<Exercise> exercise =
        boost::make_shared<BermudanExercise>(exerciseDates, false);
    boost::shared_ptr<Swaption_t<dbl> > swaption =
        boost::make_shared<Swaption_t<dbl> >(underlying, exercise);
    boost::shared_ptr<Swaption_t<dblAD> > swaptionAD =
        boost::make_shared<Swaption_t<dblAD> >(underlyingAD, exercise);

    // the gsr model

    std::vector<Date> stepDates(exerciseDates.begin(), exerciseDates.end() - 1);

    // standard sigmas
    std::vector<dbl> sigmas(stepDates.size() + 1, 0.01);

    // sigmas calibrated
    // std::vector<dbl> sigmas;
    // strike 5 percent
    // sigmas += 0.0048487, 0.00486706, 0.00489049, 0.00487795, 0.00492656,
    //     0.00490272, 0.00493519, 0.00492715, 0.0049653;
    // strike atm
    // sigmas += 0.00372794749994, 0.00374715655346, 0.00374209625552,
    //     0.0037592713444, 0.00379849214083, 0.00375681533587,
    //     0.00377164551162,
    //     0.00378151210676, 0.00385943120977;
    // pertube sigmas, so to get sensis
    // for (Size i = 0; i < sigmas.size(); ++i)
    //     sigmas[i] *= 1+1E-6;

    dbl reversion = 0.01;

    std::vector<dblAD> sigmasAD(sigmas.begin(), sigmas.end());

    // Compute Sensis to sigmas

    // CppAD::Independent(sigmasAD);

    boost::shared_ptr<Gsr_t<dbl> > gsr =
        boost::make_shared<Gsr_t<dbl> >(yts6m_h, stepDates, sigmas, reversion);
    boost::shared_ptr<Gsr_t<dblAD> > gsrAD = boost::make_shared<Gsr_t<dblAD> >(
        yts6mAD_h, stepDates, sigmasAD, reversion);

    boost::shared_ptr<PricingEngine> swaptionEngine =
        boost::make_shared<Gaussian1dSwaptionEngine_t<dbl> >(gsr,32,5.0);
    boost::shared_ptr<PricingEngine> swaptionEngineAD =
        boost::make_shared<Gaussian1dSwaptionEngine_t<dblAD> >(gsrAD,32,5.0);

    // calibration basket

    std::vector<dbl> inputVol(9, 0.0);
    std::vector<dblAD> inputVolAD(9, 0.0);

    std::vector<boost::shared_ptr<CalibrationHelper_t<dbl> > > basket;
    std::vector<boost::shared_ptr<CalibrationHelper_t<dblAD> > > basketAD;
    for (Size i = 1; i < 10; ++i) {
        Period swapLength = (10 - i) * Years;
        inputVol[i - 1] =
            swvol->volatility(exerciseDates[i - 1], swapLength, strike);
        inputVolAD[i - 1] = inputVol[i - 1];
    }

    CppAD::Independent(inputVolAD);

    for (Size i = 1; i < 10; ++i) {
        Period swapLength = (10 - i) * Years;
        Handle<Quote_t<dbl> > vol_q(
            boost::make_shared<SimpleQuote_t<dbl> >(inputVol[i - 1]));
        Handle<Quote_t<dblAD> > volAD_q(
            boost::make_shared<SimpleQuote_t<dblAD> >(inputVolAD[i - 1]));
        boost::shared_ptr<SwaptionHelper_t<dbl> > tmp =
            boost::make_shared<SwaptionHelper_t<dbl> >(
                exerciseDates[i - 1], swapLength, vol_q, euribor6m, 1 * Years,
                Thirty360(), Actual360(), yts6m_h,
                CalibrationHelper_t<dbl>::RelativePriceError, strike, 1.0);
        tmp->setPricingEngine(swaptionEngine);
        basket.push_back(tmp);
        boost::shared_ptr<SwaptionHelper_t<dblAD> > tmpAD =
            boost::make_shared<SwaptionHelper_t<dblAD> >(
                exerciseDates[i - 1], swapLength, volAD_q, euribor6mAD,
                1 * Years, Thirty360(), Actual360(), yts6mAD_h,
                CalibrationHelper_t<dblAD>::RelativePriceError, strike, 1.0);
        tmpAD->setPricingEngine(swaptionEngineAD);
        basketAD.push_back(tmpAD);
    }

    // calibrate the model

    LevenbergMarquardt_t<dbl> method;
    EndCriteria_t<dbl> ec(1000, 10, 1E-8, 1E-8, 1E-8);
    LevenbergMarquardt_t<dblAD> methodAD;
    EndCriteria_t<dblAD> ecAD(1000, 10, 1E-8, 1E-8, 1E-8);

    timer.start();
    gsr->calibrateVolatilitiesIterative(basket, method, ec);
    timer.stop();
    std::cout << "dbl model calibration timing = " << timer.elapsed()
              << std::endl;

    timer.start();
    gsrAD->calibrateVolatilitiesIterative(basketAD, methodAD, ecAD);
    timer.stop();
    std::cout << "AD model calibration timing = " << timer.elapsed()
              << std::endl;

    Array_t<dbl> volatility = gsr->volatility();
    std::cout << "Expiry;Sigma;PriceMdl;PriceMkt;MktVol" << std::endl;
    for (Size i = 0; i < basket.size(); ++i) {
        boost::shared_ptr<SwaptionHelper_t<dbl> > helper =
            boost::dynamic_pointer_cast<SwaptionHelper_t<dbl> >(basket[i]);
        std::cout << std::setprecision(12)
                  << helper->swaption()->exercise()->date(0) << ";"
                  << volatility[i] << ";" << helper->modelValue() << ";"
                  << helper->marketValue() << ";" << inputVol[i] << std::endl;
    }
    std::cout << std::endl;

    Array_t<dblAD> volatilityAD = gsrAD->volatility();
    std::cout << "Expiry;Sigma;PriceMdl;PriceMkt;MktVol" << std::endl;
    for (Size i = 0; i < basket.size(); ++i) {
        boost::shared_ptr<SwaptionHelper_t<dblAD> > helper =
            boost::dynamic_pointer_cast<SwaptionHelper_t<dblAD> >(basketAD[i]);
        std::cout << helper->swaption()->exercise()->date(0) << ";"
                  << volatilityAD[i] << ";" << helper->modelValue() << ";"
                  << helper->marketValue() << ";" << inputVolAD[i] << std::endl;
    }
    std::cout << std::endl;

    // price the calibration helpers and compute their dNPV / dsigma

    // timer.start();

    //  std::vector<dblAD> basketNpv(basketAD.size(), 0.0);

    //  for (Size i = 0; i < basketAD.size(); ++i) {
    //      basketAD[i]->setPricingEngine(swaptionEngineAD);
    //      basketNpv[i] = basketAD[i]->modelValue();
    //      std::cout << "basket #" << i << " npv = " << basketNpv[i] <<
    //      std::endl;
    //  }

    //  CppAD::ADFun<Real> f(sigmasAD, basketNpv);
    //  std::vector<Real> helperVega(basketAD.size() * sigmasAD.size());
    //  helperVega = f.Jacobian(sigmas);
    //  Matrix_t<dbl> helperVega_m(basketAD.size(), sigmasAD.size(),
    //                               helperVega.begin(), helperVega.end());

    //  std::cout << "Jacobian of swaption helpers by sigma" << std::endl;
    //  for(Size i=0;i<basketAD.size();++i) {
    //      for(Size j=0;j<sigmasAD.size();++j) {
    //          std::cout << helperVega_m[i][j] << ";";
    //      }
    //      std::cout << std::endl;
    //  }

    // // output the helpers vega
    //  for (Size i = 0; i < basket.size(); ++i) {
    //      boost::shared_ptr<Swaption_t<dbl>> sw =
    //      boost::dynamic_pointer_cast<SwaptionHelper_t<dbl>>(basket[i])->swaption();
    //      boost::shared_ptr<BlackSwaptionEngine_t<dbl>> en =
    //      boost::make_shared<BlackSwaptionEngine_t<dbl>>(yts6m_h,
    //      inputVol[i]);
    //      sw->setPricingEngine(en);
    //      Real vega = sw->result<Real>("vega");
    //      std::cout << "helper vega #" << i << " = " << vega << std::endl;
    //  }

    //  timer.stop();
    //  std::cout << "additional stuff = " << timer.elapsed() << std::endl;

    // price the bermudan and sensis to sigma

    timer.start();
    swaption->setPricingEngine(swaptionEngine);
    std::cout << "berm npv=" << swaption->NPV() << std::endl;
    timer.stop();
    std::cout << "plain pricing = " << timer.elapsed() << std::endl;

    timer.start();
    swaptionAD->setPricingEngine(swaptionEngineAD);
    std::vector<dblAD> yAD(1, 0.0);
    yAD[0] = swaptionAD->NPV();
    std::cout << "Bermudan NPV = " << std::setprecision(12) << yAD[0]
              << std::endl;
    timer.stop();
    std::cout << "AD pricing =" << timer.elapsed() << std::endl;

    timer.start();
    CppAD::ADFun<Real> f(inputVolAD, yAD);
    // CppAD::ADFun<Real> f(sigmasAD, yAD);
    std::vector<Real> vega(sigmasAD.size()), w(1, 1.0);
    vega = f.Reverse(1, w);
    timer.stop();
    std::cout << "AD Deltas: " << timer.elapsed() << std::endl;

    for (Size i = 0; i < vega.size(); ++i)
        std::cout << "vega #" << i << " = " << vega[i] << std::endl;

    return 0;
}
