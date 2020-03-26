package main

import (
	"encoding/json"
	"fmt"
	"log"
	"net/http"
)

// Experiment defines an experiment structure
type Experiment struct {
	Name string `json:"name"`
	Date string `json:"date"`
	// key value parameters
	Parmeters   map[string]interface{} `json:"parameters"`
	Successfull bool                   `json:"successful"`
}

// ExperimentSlice holds lists of experiments
type ExperimentSlice struct {
	Experiments     []Experiment `json:"experiments"`
	ExperimentCount int          `json:"experimentCount"`
}

// experiments holds the list of all experiments
var experiments ExperimentSlice

func getParticularExperiment(queryType, queyrValue string) (*Experiment, error) {
	switch queryType {
	case "NAME":
		for _, exp := range experiments.Experiments {
			if exp.Name == queryType {
				return &exp, nil
			}
		}
	case "DATE":
		for _, exp := range experiments.Experiments {
			if exp.Date == queryType {
				return &exp, nil
			}
		}
	}

	return nil, fmt.Errorf("Not found %s in query: %s ", queyrValue, queryType)
}

func handleExperimentSubmission(w http.ResponseWriter, r *http.Request) {

	switch r.Method {
	case "GET":
		for k, v := range r.URL.Query() {
			fmt.Printf("%s: %s\n", k, v)
		}
		w.Write([]byte("Received a GET request\n"))
	case "POST":
		// retrieve the message body
		var exp Experiment

		err := json.NewDecoder(r.Body).Decode(&exp)
		if err != nil {
			http.Error(w, err.Error(), http.StatusBadRequest)
			return
		}
		experiments.Experiments = append(experiments.Experiments, exp)
		experiments.ExperimentCount++
		w.Write([]byte("Added an experiment to the experiment list\n"))
	default:
		w.WriteHeader(http.StatusNotImplemented)
		w.Write([]byte(http.StatusText(http.StatusNotImplemented)))
	}

}

func experimentsToJSON() ([]byte, error) {
	experimentJSON, err := json.Marshal(experiments)
	return experimentJSON, err
}

func handleListExperiments(w http.ResponseWriter, r *http.Request) {
	switch r.Method {
	case "GET":

		expJ, err := experimentsToJSON()
		if err != nil {
			log.Panic(err)
		}

		w.Header().Set("Content-Type", "application/json")
		w.WriteHeader(http.StatusOK)

		w.Write(expJ)
	}
}

func main() {

	experiments = ExperimentSlice{
		ExperimentCount: 0,
	}

	mux := http.NewServeMux()
	mux.HandleFunc("/", handleExperimentSubmission)
	mux.HandleFunc("/experiments", handleListExperiments)

	log.Println("Starting server on :4000...")
	err := http.ListenAndServe(":4000", mux)
	log.Fatal(err)
}
