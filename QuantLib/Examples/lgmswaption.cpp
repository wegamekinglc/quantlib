#include <ql/quantlib.hpp>

using namespace QuantLib;

int main() {

    try {

        Real rateLevel = 0.02;
        Date evalDate(12, January, 2015);

        Settings::instance().evaluationDate() = evalDate;

        std::clog << "reference date " << evalDate << std::endl;

        Handle<YieldTermStructure> yts(boost::make_shared<FlatForward>(
            evalDate, rateLevel, Actual365Fixed()));

        boost::shared_ptr<IborIndex> euribor6m =
            boost::make_shared<Euribor>(6 * Months, yts);

        std::vector<boost::shared_ptr<Swaption> > swaptions;

        std::vector<Date> exerciseDates;

        for (Size i = 0; i < 10; ++i) {

            Real strike = 0.02 + static_cast<Real>(i) * 0.0010;

            Date effectiveDate = TARGET().advance(evalDate, 2 * Days);
            Date startDate = TARGET().advance(effectiveDate, 1 * Years);
            Date maturityDate =
                TARGET().advance(startDate, 10 * Years); // maturity here !!!

            Schedule fixedSchedule(startDate, maturityDate, 1 * Years, TARGET(),
                                   ModifiedFollowing, ModifiedFollowing,
                                   DateGeneration::Forward, false);
            Schedule floatingSchedule(startDate, maturityDate, 6 * Months,
                                      TARGET(), ModifiedFollowing,
                                      ModifiedFollowing,
                                      DateGeneration::Forward, false);

            // std::clog << "exercise dates:" << std::endl;
            exerciseDates.clear();
            for (Size ii = 0; ii < 9; ++ii) {
                exerciseDates.push_back(
                    TARGET().advance(fixedSchedule[ii], -2 * Days));
                // std::clog << exerciseDates.back() << "\n";
            }

            boost::shared_ptr<Exercise> exercise =
                boost::make_shared<BermudanExercise>(exerciseDates, false);

            boost::shared_ptr<VanillaSwap> underlying =
                boost::make_shared<VanillaSwap>(VanillaSwap(
                    VanillaSwap::Payer, 1.0, fixedSchedule, strike, Thirty360(),
                    floatingSchedule, euribor6m, 0.0, Actual360()));

            boost::shared_ptr<Swaption> swaption =
                boost::make_shared<Swaption>(underlying, exercise);

            swaptions.push_back(swaption);
        }

        std::vector<Date> stepDates(exerciseDates.begin(),
                                    exerciseDates.end() - 1);

        std::vector<Real> sigmas(stepDates.size() + 1, 0.0050);
        Real reversion = 0.0;

        boost::shared_ptr<Gsr> gsr =
            boost::make_shared<Gsr>(yts, stepDates, sigmas, reversion, 50.0);

        boost::shared_ptr<Lgm1> lgm =
            boost::make_shared<Lgm1>(yts, stepDates, sigmas, reversion);

        boost::shared_ptr<PricingEngine> swaptionEngineGsr =
            boost::make_shared<Gaussian1dSwaptionEngine>(gsr, 64, 6.0, false,
                                                         false);

        boost::shared_ptr<PricingEngine> swaptionEngineLgm =
            boost::make_shared<Gaussian1dSwaptionEngine>(lgm, 64, 6.0, false,
                                                         false);

        // ----------------------------------
        // test fortran (ad) engine
        // ----------------------------------

        boost::shared_ptr<PricingEngine> swaptionEngineLgmAD =
            boost::make_shared<LgmSwaptionEngineAD>(lgm, 64, 6.0);

        for (Size i = 0; i < 1; ++i) {
            swaptions[i]->setPricingEngine(swaptionEngineLgmAD);
            Real npvLgmAD = swaptions[i]->NPV();
            std::clog << "i=" << i << " npv (LgmAD) = " << npvLgmAD
                      << std::endl;
            std::vector<Real> times =
                swaptions[i]->result<std::vector<Real> >("sensitivityTimes");
            std::vector<Real> sensH =
                swaptions[i]->result<std::vector<Real> >("sensitivityH");
            std::vector<Real> sensZ =
                swaptions[i]->result<std::vector<Real> >("sensitivityZeta");
            std::vector<Real> sensD =
                swaptions[i]->result<std::vector<Real> >("sensitivityDiscount");
            for (Size j = 0; j < times.size(); ++j) {
                std::clog << j << ";" << times[j] << ";" << sensH[j] << ";"
                          << sensZ[j] << ";" << sensD[j] << std::endl;
            }
        }

        return 0;

    } catch (QuantLib::Error e) {
        std::clog << e.what() << std::endl;
    } catch (std::exception e) {
        std::clog << e.what() << std::endl;
    }

} // main
