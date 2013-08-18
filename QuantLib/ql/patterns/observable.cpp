/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2013 Klaus Spanderen

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

/*! \file observable.cpp
    \brief observer/observable pattern
*/

#include <ql/patterns/observable.hpp>
#include <boost/thread/lock_guard.hpp>

namespace QuantLib {
    class Observer::Proxy {
      public:
        Proxy(Observer* observer)
        : active_(true),
          observer_(observer) {
        }
        void update() const {
            boost::lock_guard<boost::mutex> lock(mutex_);
            if(active_) {
                observer_->update();
            }
        }
        void deactivate() {
            boost::lock_guard<boost::mutex> lock(mutex_);
            active_ = false;
        }

        bool isActive() const {
            return active_;
        }
      private:
        bool active_;
        mutable boost::mutex mutex_;
        Observer* const observer_;
    };

    void Observer::deactivate() {
        proxy_->deactivate();
    }

    bool Observer::isActive() const {
        return proxy_->isActive();
    }

    Observer::Observer()
    : proxy_(new Proxy(this)) {
    }

    std::pair<std::set<boost::shared_ptr<Observable> >::iterator, bool>
    Observer::registerWith(const boost::shared_ptr<Observable>& h) {
        boost::lock_guard<boost::mutex> lock(mutex_);
        if (h) {
            h->registerObserver(proxy_);
            return observables_.insert(h);
        }
        return std::make_pair(observables_.end(), false);
    }


    Size Observer::unregisterWith(const boost::shared_ptr<Observable>& h) {
        boost::lock_guard<boost::mutex> lock(mutex_);
        if (h) {
            h->unregisterObserver(proxy_);
        }

        return observables_.erase(h);
    }

    Observer::Observer(const Observer& o)
    : proxy_(new Proxy(this)) {
        {
            boost::lock_guard<boost::mutex> lock(o.mutex_);
            observables_ = o.observables_;
        }
        for (iterator i=observables_.begin(); i!=observables_.end(); ++i)
            (*i)->registerObserver(proxy_);
    }

    Observer& Observer::operator=(const Observer& o) {
        boost::lock_guard<boost::mutex> lock(mutex_);
        for (iterator i=observables_.begin(); i!=observables_.end(); ++i)
            (*i)->unregisterObserver(proxy_);

        {
            boost::lock_guard<boost::mutex> lock(o.mutex_);
            observables_ = o.observables_;
        }
        for (iterator i=observables_.begin(); i!=observables_.end(); ++i)
            (*i)->registerObserver(proxy_);
        return *this;
    }

    Observer::~Observer() {
        for (iterator i=observables_.begin(); i!=observables_.end(); ++i)
            (*i)->unregisterObserver(proxy_);
    }


    void Observable::registerObserver(
        const boost::shared_ptr<Observer::Proxy>& observerProxy) {

        signal_type::slot_type slot(&Observer::Proxy::update,
                                    observerProxy.get());
        sig_.connect(slot.track(observerProxy));
    }

    void Observable::unregisterObserver(
        const boost::shared_ptr<Observer::Proxy>& observerProxy) {

        sig_.disconnect(boost::bind(&Observer::Proxy::update,
                        observerProxy.get()));
    }

    void Observable::notifyObservers() {
        sig_();
    }
}
