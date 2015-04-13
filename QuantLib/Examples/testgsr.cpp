#include <ql/quantlib.hpp>

#include <boost/assign/std/vector.hpp>
#include <boost/make_shared.hpp>
#include <boost/timer.hpp>

using namespace boost::assign;
using namespace QuantLib;

// typedef Real dbl;
typedef CppAD::AD<double> dbl;

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
    boost::shared_ptr<Euribor_t<dbl> > euribor6mBt =
        boost::make_shared<Euribor_t<dbl> >(6 * Months);
    std::vector<boost::shared_ptr<RateHelper_t<dbl>::Type> > helpers;
    for (Size i = 1; i <= 30; ++i) {
        boost::shared_ptr<RateHelper_t<dbl>::Type> tmp =
            boost::make_shared<SwapRateHelper_t<dbl> >(
                swapQuote, i * Years, TARGET(), Annual, ModifiedFollowing,
                Thirty360(), euribor6mBt);
        helpers.push_back(tmp);
    }
    boost::shared_ptr<YieldTermStructure_t<dbl> > yts6m = boost::make_shared<
        PiecewiseYieldCurve<ZeroYield, Linear, IterativeBootstrap, dbl> >(
        refDate, helpers, Actual365Fixed());

    Handle<YieldTermStructure_t<dbl> > yts6m_h(yts6m);

    yts6m_h->enableExtrapolation();

    std::cout << "test yts = " << yts6m_h->discount(1.0) << std::endl;
    timer.stop();
    std::cout << "timing: " << timer.elapsed() << std::endl;

    // index

    timer.start();

    boost::shared_ptr<Euribor_t<dbl> > euribor6m =
        boost::make_shared<Euribor_t<dbl> >(6 * Months, yts6m_h);

    boost::shared_ptr<Euribor_t<dbl> > euribor6mTmp =
        boost::make_shared<Euribor_t<dbl> >(6 * Months, yts6m_h);

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

    dbl h = 1E-4; // for fd approximation
    std::vector<dbl> xAD;
    xAD += dbl(0.03);
    xAD += dbl(0.60);
    xAD += dbl(0.12);
    xAD += dbl(0.30);

    CppAD::Independent(xAD);

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
            boost::make_shared<SimpleQuote_t<dbl> >(xAD[0])); // alpha
        sabrTmp += Handle<Quote_t<dbl> >(
            boost::make_shared<SimpleQuote_t<dbl> >(xAD[1])); // beta
        sabrTmp += Handle<Quote_t<dbl> >(
            boost::make_shared<SimpleQuote_t<dbl> >(xAD[2])); // nu
        sabrTmp += Handle<Quote_t<dbl> >(
            boost::make_shared<SimpleQuote_t<dbl> >(xAD[3])); // rho
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

    dbl strike = 0.05;

    boost::shared_ptr<VanillaSwap_t<dbl> > underlying =
        boost::make_shared<VanillaSwap_t<dbl> >(
            VanillaSwap_t<dbl>::Payer, 1.0, fixedSchedule, strike, Thirty360(),
            floatingSchedule, euribor6m, 0.0, Actual360());

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

    // the gsr model

    std::vector<Date> stepDates(exerciseDates.begin(), exerciseDates.end() - 1);

    // standard sigmas
    //std::vector<dbl> sigmas(stepDates.size() + 1, 0.0050);

    // sigmas calibrated
    std::vector<dbl> sigmas;
    sigmas += 0.0048487, 0.00486706, 0.00489049, 0.00487795, 0.00492656,
        0.00490272, 0.00493519, 0.00492715, 0.0049653;
    // pertube sigmas, so to get sensis
    for(Size i=0;i<sigmas.size();++i)
        sigmas *= 1.0001;

    dbl reversion = 0.01;

    boost::shared_ptr<Gsr_t<dbl> > gsr =
        boost::make_shared<Gsr_t<dbl> >(yts6m_h, stepDates, sigmas, reversion);

    boost::shared_ptr<PricingEngine> swaptionEngine =
        boost::make_shared<Gaussian1dSwaptionEngine_t<dbl> >(gsr);

    swaption->setPricingEngine(swaptionEngine);

    // calibration basket

    std::vector<boost::shared_ptr<CalibrationHelper_t<dbl> > > basket;
    for (Size i = 1; i < 10; ++i) {
        Period swapLength = (10 - i) * Years;
        dbl vol = swvol->volatility(exerciseDates[i - 1], swapLength, strike);
        Handle<Quote_t<dbl> > vol_q(
            boost::make_shared<SimpleQuote_t<dbl> >(vol));
        boost::shared_ptr<SwaptionHelper_t<dbl> > tmp =
            boost::make_shared<SwaptionHelper_t<dbl> >(
                exerciseDates[i - 1], swapLength, vol_q, euribor6m, 1 * Years,
                Thirty360(), Actual360(), yts6m_h,
                CalibrationHelper_t<dbl>::RelativePriceError, strike, 1.0);
        tmp->setPricingEngine(swaptionEngine);
        basket.push_back(tmp);
    }

    // calibrate the model

    LevenbergMarquardt_t<dbl> method;
    EndCriteria_t<dbl> ec(1000, 10, 1E-8, 1E-8, 1E-8);

    gsr->calibrateVolatilitiesIterative(basket, method, ec);

    Array_t<dbl> volatility = gsr->volatility();

    std::cout << "Expiry;Sigma;PriceMdl;PriceMkt" << std::endl;
    for (Size i = 0; i < basket.size(); ++i) {
        boost::shared_ptr<SwaptionHelper_t<dbl> > helper =
            boost::dynamic_pointer_cast<SwaptionHelper_t<dbl> >(basket[i]);
        std::cout << helper->swaption()->exercise()->date(0) << ";"
                  << volatility[i] << ";" << helper->modelValue() << ";"
                  << helper->marketValue() << std::endl;
    }
    std::cout << std::endl;

    timer.stop();
    std::cout << "timing: " << timer.elapsed() << std::endl;

    // price the bermudan

    timer.start();

    std::vector<dbl> yAD(1, 0.0);

    yAD[0] = swaption->NPV();

    CppAD::ADFun<Real> f(xAD, yAD);
    std::vector<Real> sabrDelta(xAD.size()), w(1, 1.0);
    sabrDelta = f.Reverse(1, w);
    timer.stop();

    std::cout << "Bermudan NPV = " << yAD[0] << std::endl;
    std::cout << "dNPV/dalpha = " << sabrDelta[0] << std::endl;
    std::cout << "dNPV/dbeta = " << sabrDelta[1] << std::endl;
    std::cout << "dNPV/dnu = " << sabrDelta[2] << std::endl;
    std::cout << "dNPV/drho = " << sabrDelta[3] << std::endl;

    std::cout << "timing: " << timer.elapsed() << std::endl;

    return 0;
}
