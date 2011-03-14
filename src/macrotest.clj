(defmacro test0 [x]
  (keyword x))

(defmacro test [lst]
  (apply hash-map
         (apply concat (for [e lst]
                         [(keyword e) (str "VK_" e)]))))
