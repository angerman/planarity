(ns utils)

(defn sleep
  ([n] (. Thread (sleep n)))
  ([] (sleep 2000)))
