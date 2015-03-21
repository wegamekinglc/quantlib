(defun templateize (classname)
"some standard query-replace commands to support teamplate'ization"
(set-buffer (current-buffer))
(interactive "sClass name: \n")
; replace standard types
(query-replace-regexp "\\(double\\|Real\\|Spread\\|Volatility\\|Rate\\)\\([ >]\\)" 
                      "T\\2" nil (point-min) (point-max))
; replace class name
(query-replace-regexp (concat classname "\\([( ]\\)") (concat classname "_t\\1") nil (point-min) (point-max))
(query-replace-regexp (concat classname "::") (concat classname "_t<T>::") nil (point-min) (point-max))
; replace frequent other class names
(query-replace-regexp "\\(Interpolation\\|Quote\\|SimpleQuote\\|SmileSection\\)\\([ (]\\)" "\\1_t<T>\\2" 
                      nil (point-min) (point-max))
; replace std:: functions
(query-replace-regexp "\\(std::\\)\*\\(max\\|min\\|pow\\|log\\|exp\\|abs\\|sqrt\\)(" 
                      "QLFCT::\\2(" nil (point-min) (point-max))
)
