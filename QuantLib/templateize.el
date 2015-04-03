(defun templateize (classname)
"some standard query-replace commands to support template'ization"
(interactive "sClass name: \n")

(set-buffer (current-buffer))

; replace standard types
(query-replace-regexp "\\([ ,<(\n\r]\\)\\(double\\|Real\\|Spread\\|Volatility\\|Rate\\)\\([ >),\n\r]\\)"
                      "\\1T\\3" nil (point-min) (point-max))

; replace class name
(query-replace-regexp (concat classname " \\([a-zA-Z]\\)")
                      (concat classname "_t<T> \\1") nil (point-min) (point-max))
(query-replace-regexp (concat classname "::")
                      (concat classname "_t<T>::") nil (point-min) (point-max))
(query-replace-regexp (concat classname "\\([( ]\\)") (concat classname "_t\\1") nil (point-min) (point-max))

; replace frequent other class names
(query-replace-regexp "\\([ ,<\n\r]\\)\\(Array\\|BilinearInterpolation\\|BackwardFlatLinearInterpolation\\|Interpolation\\|Interpolation2D\\|Matrix\\|Quote\\|SABRInterpolation\\|SimpleQuote\\|SmileSection\\|SwapIndex\\|VanillaSwap\\|VolatilityTermStructure\\)\\([ (>\n\r\]\\)" "\\1\\2_t<T>\\3"
                      nil (point-min) (point-max))

; replace std:: functions
(query-replace-regexp "\\(std::\\)\*fabs(" "QLFCT::abs(")
(query-replace-regexp "\\(std::\\)\*\\(abs\\|max\\|min\\|pow\\|log\\|exp\\|sqrt\\)("
                      "QLFCT::\\2(" nil (point-min) (point-max))
)
