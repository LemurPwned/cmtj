package main

import (
	"errors"
	"log"
)

func findExperimentByID(experiments *ExperimentSlice, experimentID string) (*Experiment, int, error) {
	log.Println("Seaching for the id", experimentID)
	for index, experiment := range experiments.Experiments {
		if experiment.ID == experimentID {
			return &experiment, index, nil
		}
	}
	return nil, -1, errors.New("No experiment found with id")
}

func addMetricsToExperiment(experiments *ExperimentSlice,
	experimentID, metricName string, metricValue float32) {
	exp, _, err := findExperimentByID(experiments, experimentID)
	if err == nil {
		log.Println(err)
	}
	exp.Metrics[metricName] = append(exp.Metrics[metricName], metricValue)
}

func deleteExperiment(experiments *ExperimentSlice, experimentID string) error {
	_, index, err := findExperimentByID(experiments, experimentID)
	if err != nil {
		log.Println(err)
		return err
	}
	l := experiments.ExperimentCount
	experiments.Experiments[index] = experiments.Experiments[l-1]
	// We do not need to put experiment[i] at the end, as it will be discarded anyway
	experiments.Experiments = experiments.Experiments[:l-1]
	experiments.ExperimentCount = l - 1
	return nil
}
