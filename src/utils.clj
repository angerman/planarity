(ns utils)

(defn none? [pred coll]
  (every? (complement pred) coll))

(defn sleep
  ([n] (. Thread (sleep n)))
  ([] (sleep 2000)))
