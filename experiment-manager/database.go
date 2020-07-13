package main

import (
	"context"
	"errors"
	"fmt"
	"log"
	"time"

	"go.mongodb.org/mongo-driver/mongo"
	"go.mongodb.org/mongo-driver/mongo/options"
	"go.mongodb.org/mongo-driver/mongo/readpref"
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

func connectMongo() {
	ctx, cancel := context.WithTimeout(context.Background(), 10*time.Second)
	defer cancel()
	client, err := mongo.Connect(ctx, options.Client().ApplyURI("mongodb://localhost:27017"))
	if err != nil {
		log.Panicln(err)
	}
	defer func() {
		if err = client.Disconnect(ctx); err != nil {
			panic(err)
		}
	}()

	if err := client.Ping(ctx, readpref.Primary()); err != nil {
		panic(err)
	}

	fmt.Println("Successfully connected and pinged.")
}
