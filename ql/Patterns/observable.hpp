
/*
 Copyright (C) 2000, 2001, 2002 RiskMap srl

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software developed by the QuantLib Group; you can
 redistribute it and/or modify it under the terms of the QuantLib License;
 either version 1.0, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 QuantLib License for more details.

 You should have received a copy of the QuantLib License along with this
 program; if not, please email ferdinando@ametrano.net

 The QuantLib License is also available at http://quantlib.org/license.html
 The members of the QuantLib Group are listed in the QuantLib License
*/
/*! \file observable.hpp
    \brief observer/observable pattern

    \fullpath
    ql/Patterns/%observable.hpp
*/

// $Id$

#ifndef quantlib_observable_h
#define quantlib_observable_h

#include <ql/qldefines.hpp>
#include <set>

namespace QuantLib {

    //! Implementations of design patterns
    /*! See sect. \ref patterns */
    namespace Patterns {

        //! Object that gets notified when a given observable changes
        /*! These classes are a simplified version of the ones implemented in
            Bruce Eckel, Thinking in C++ (http://www.bruceeckel.com) which in
            turn mirror the Java Observer and Observable interfaces.
        */
        class Observer {
          public:
            virtual ~Observer() {}
            /*! This method must be implemented in derived classes. An
                instance of %Observer does not call this method directly:
                instead, it will be called by the observables the instance
                registered with when they need to notify any changes.
            */
            virtual void update() = 0;
        };

        //! Object that notifies its changes to a set of observables
        class Observable {
          public:
            virtual ~Observable() {}
            /*! \name Observer management
                \warning It is responsibility of the programmer to make sure
                that the registered observers unregister themselves before
                going out of scope.
            */
            //@{
            virtual void registerObserver(Observer*);
            void registerObservers(std::set<Observer*>&);
            virtual void unregisterObserver(Observer*);
            void unregisterObservers(std::set<Observer*>&);
            virtual void unregisterAll();
            std::set<Observer*> observers() const;
            /*! In derived classes, this method should be called at the end
                of non-const methods or when the programmer desires to notify
                any changes.
            */
            virtual void notifyObservers();
            //@}
          private:
            std::set<Observer*> theObservers;
        };


        // inline definitions

        inline void Observable::registerObserver(Observer* o) {
            theObservers.insert(o);
        }

        inline void Observable::registerObservers(std::set<Observer*>& o) {
            for (std::set<Observer*>::iterator i = o.begin(); i != o.end(); ++i)
                registerObserver(*i);
        }

        inline void Observable::unregisterObserver(Observer* o) {
            theObservers.erase(o);
        }

        inline void Observable::unregisterObservers(std::set<Observer*>& o) {
            for (std::set<Observer*>::iterator i = o.begin(); i != o.end(); ++i)
                unregisterObserver(*i);
        }

        inline void Observable::unregisterAll() {
            theObservers.clear();
        }

        inline std::set<Observer*> Observable::observers() const {
            return theObservers;
        }

        inline void Observable::notifyObservers() {
            for (std::set<Observer*>::iterator i = theObservers.begin();
                i != theObservers.end(); ++i)
                    (*i)->update();
        }

    }

}


#endif
