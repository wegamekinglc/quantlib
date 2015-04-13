(defun templateize (classname)
"some standard query-replace commands to support template'ization"
(interactive "sClass name: \n")

(set-buffer (current-buffer))

; replace standard types
(query-replace-regexp "\\([ ,<(\n\r]\\)\\(double\\|Real\\|Spread\\|Volatility\\|Rate\\)\\([ *>),\n\r]\\)"
                      "\\1T\\3" nil (point-min) (point-max))

; replace class name
(query-replace-regexp (concat classname "\\( {\\|(\\)") (concat classname "_t\\1") nil (point-min) (point-max))
(query-replace-regexp (concat classname "\\( \\|::\\|>\\)")
                      (concat classname "_t<T>\\1") nil (point-min) (point-max))

; replace frequent other class names
(query-replace-regexp "\\([ (,<\n\r]\\)\\(Array\\|BilinearInterpolation\\|BackwardFlatLinearInterpolation\\|BlackCapFloorEngine\\|BlackSwaptionEngine\\|CalibrationFunction\\|CalibrationHelper\\|CapFloor\\|Constraint\\|CostFunction\\|CumulativeNormalDistribution\\|DiscretizedAsset\\|EndCriteria\\|FixedRateCoupon\\|FloatingRateCoupon\\|Gaussian1dModel\\|IborIndex\\|Instrument\\|InterestRate\\|InterestRateIndex\\|Interpolation\\|Interpolation2D\\|Leg\\|Matrix\\|NewtonSafe\\|OptimizationMethod\\|Problem\\|Projection\\|ProjectedConstraint\\|ProjectedCostFunction\\|Quote\\|SABRInterpolation\\|SimpleQuote\\|SmileSection\\|SwapIndex\\|Swaption\\|SwaptionVolatilityStructure\\|TridiagonalOperator\\|VanillaSwap\\|VolatilityTermStructure\\|YieldTermStructure\\)\\([ :&(>\n\r\]\\)" "\\1\\2_t<T>\\3"
                      nil (point-min) (point-max))

; add reference for frequent base class members and methods
(query-replace-regexp "\\([ *([\r\n]\\)\\(calculate\\|registerWith\\|notifiyObservers\\|arguments_\\|results_\\)" "\\1this->\\2" nil (point-min) (point-max))

; replace std:: functions
(query-replace-regexp "\\(std::\\)\*fabs(" "QLFCT::abs(")
(query-replace-regexp "\\(std::\\)\*\\(abs\\|max\\|min\\|pow\\|log\\|exp\\|sqrt\\|sin\\|cos\\|tan\\|asin\\|acos\\|atan\\|sinh\\|cosh\\\tanh\\)("
                      "QLFCT::\\2(" nil (point-min) (point-max))

)
