/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2013 Jose Aparicio

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#ifndef quantlib_credit_index_hpp
#define quantlib_credit_index_hpp

#include <ql/index.hpp>
#include <ql/time/calendar.hpp>
#include <ql/currency.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/time/period.hpp>
#include <ql/handle.hpp>
#include <ql/experimental/credit/defaulttype.hpp>
#include <ql/experimental/credit/issuer.hpp>
#include <ql/experimental/credit/recoveryratequote.hpp>
#include <ql/time/dategenerationrule.hpp>

namespace QuantLib {

    class DefaultProbabilityTermStructure;
    class YieldTermStructure;
    class CreditDefaultSwap;
    /*
    \to do Need a more elaborated fixing date mechanism
    For example (see markit web page documentation) iTraxx "Fixings are \
    set at 11am London time every Friday, and at 4pm on every IMM roll date: \
    20th March; June; September and December."

    Maybe nothing should be set at this level and leave the fixings dates \
    to derived classes.

    See example: http://www.creditfixings.com/CreditEventAuctions/fixings.jsp

    \todo
    Yet one can argue that an index fixing is a quote of two numbers; the spread
    itself and the recovery. As it is today IndexManager can store Reals only
    and not an arbitrary type or agregation types. I think this is the base of 
    some potential inconsistent features of the way credit coupons are priced.
    */
	class CreditIndex : public Index, public Observer {
    public:
        CreditIndex(const std::string& familyName,
                    const Period& tenor,
                    Natural fixingDays,
                    DateGeneration::Rule dateRule,
                    const Currency& currency,
                    const Calendar& fixingCalendar,
                    const DayCounter& dayCounter,
                    Seniority indexSeniority = NoSeniority);
        //! \name Index interface
        //@{
        std::string name() const;
        Calendar fixingCalendar() const;
        virtual bool isValidFixingDate(const Date& fixingDate) const = 0;
        virtual Rate fixing(const Date& fixingDate,
                            bool forecastTodaysFixing = false) const = 0;
        void addFixing(const Date& fixingDate,
                   Real fixing,
                   bool forceOverwrite);
        //@}
        virtual Handle<YieldTermStructure> nominalTermStructure() const = 0;
        //! \name Inspectors
        //@{
        std::string familyName() const;
        Period tenor() const;
        Natural fixingDays() const;
        Date fixingDate(const Date& valueDate) const;
        Seniority seniority() const {return indexSeniority_;}
        const DayCounter& dayCounter() const;
        DateGeneration::Rule rule() const;
        //@}
        /*! \name Date calculations
            @{
        */
        virtual Date valueDate(const Date& fixingDate) const;
        virtual Date maturityDate(const Date& valueDate) const;
        // @}
    protected:
        virtual Rate forecastFixing(const Date& fixingDate) const = 0;
        std::string familyName_;
        Seniority indexSeniority_;
        Period tenor_;
        Natural fixingDays_;
        Calendar fixingCalendar_;
        DayCounter dayCounter_;
        DateGeneration::Rule dateRule_;
    };

    // inline definitions

    inline std::string CreditIndex::familyName() const {
        return familyName_;
    }

    inline std::string CreditIndex::name() const {
        std::ostringstream out;
        out << tenor_;
        return familyName_ + out.str();    
    }

    inline Period CreditIndex::tenor() const {
        return tenor_;
    }

    inline Natural CreditIndex::fixingDays() const {
        return fixingDays_;
    }

    inline const DayCounter& CreditIndex::dayCounter() const {
        return dayCounter_;
    }

    inline Calendar CreditIndex::fixingCalendar() const {
        return fixingCalendar_;
    }

    inline DateGeneration::Rule CreditIndex::rule() const {
        return dateRule_;
    }


    //! Credit Default Swap index for a single underlying reference name.
    /* \todo use optional variables to specify an underlying with both \
            running and upfront
    */
    class SingleNameCreditIndex : public CreditIndex {
    public:
        SingleNameCreditIndex(
            const std::string& familyName,
            // turn to a Handle if/when Issuer becomes observable
            const Issuer& issuer, 
            const DefaultProbKey& defKey, // standard
            const Period& tenor, 
            const Frequency payFreq,  // Quarterly
            Natural fixingDays, // 0 
            DateGeneration::Rule dateRule, // CDS
            const Calendar& fixingCalendar,
            const DayCounter& dayCounter,
            const Handle<YieldTermStructure>& yts/*,
            Seniority indexSeniority*/
        ) : CreditIndex(familyName, tenor, fixingDays, 
                        dateRule, defKey.currency(), 
                        fixingCalendar, dayCounter, 
                        defKey.seniority() ), /////////indexSeniority), 
            defKey_(defKey),
            creditCurve_(issuer.defaultProbability(defKey)),
            discountTermStructure_(yts),
            issuer_(issuer),
            paymentFreq_(payFreq)
        {
            registerWith(creditCurve_);
            registerWith(discountTermStructure_);
        }

        void update();
        const DefaultProbKey& defaultKey() const;
        // defaults are to be accessed through the issuer
        const Issuer& issuer() const;
        //! \name Index interface
        //@{
        bool isValidFixingDate(const Date& fixingDate) const;
        Rate fixing(const Date& fixingDate,
                    bool forecastTodaysFixing = false) const;
        //@}
        Handle<DefaultProbabilityTermStructure> 
            defaultProbTermStructure() const;

        Handle<YieldTermStructure> nominalTermStructure() const;
        /*! Returns the zero upfront version of the underlying swap.

            \warning Relinking the term structure or the default curve 
                     underlying the index will not have effect on the 
                     returned swap.
        */
        boost::shared_ptr<CreditDefaultSwap> underlyingSwap(
                                                const Date& fixingDate) const;
    protected:
        //! Forecasts the fixing according to the relevant TS using the 
        //  conventional recovery rate
        virtual Rate forecastFixing(const Date& fixingDate) const;
        //! Forecasts the fixing according to the relevant TS
        virtual Rate forecastFixing(const Date& fixingDate, 
            Real recovery) const;
        DefaultProbKey defKey_;
        Handle<DefaultProbabilityTermStructure> creditCurve_;
        Handle<YieldTermStructure> discountTermStructure_;
        Issuer issuer_;
        Frequency paymentFreq_;
    };

    // inlines
    inline void SingleNameCreditIndex::update() {
        notifyObservers();
    }

    inline bool SingleNameCreditIndex::isValidFixingDate(
        const Date& ) const {
        return true;
    }

    inline const DefaultProbKey& SingleNameCreditIndex::defaultKey() const {
        return defKey_;
    }

    inline const Issuer& SingleNameCreditIndex::issuer() const {
        return issuer_;
    }

    inline Handle<DefaultProbabilityTermStructure> 
        SingleNameCreditIndex::defaultProbTermStructure() const {
        return creditCurve_;
    }

    inline Handle<YieldTermStructure> 
        SingleNameCreditIndex::nominalTermStructure() const {
        return discountTermStructure_;
    }

    inline Rate SingleNameCreditIndex::forecastFixing(const Date& fixingDate
        ) const {
        return forecastFixing(fixingDate, 
            RecoveryRateQuote::conventionalRecovery(indexSeniority_));
    }

}

#endif
