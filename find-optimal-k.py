#!/usr/bin/env python3

import matplotlib as mpl
mpl.use('Agg')

import argparse
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

def handler():
    parser = argparse.ArgumentParser(description="Find optimal K using elbow and silhouette methods.")
    parser.add_argument("-i", "--input_file", type=str, required=True, help="CSV file containing PC1 and PC2 columns.")
    parser.add_argument("-t", "--tag", type=str, required=True, help="Tag name for output files.")
    parser.add_argument("-kmax", "--max_clusters", type=int, default=10, help="Maximum number of clusters to test.")
    return parser.parse_args()


if __name__ == "__main__":
    args = handler()

    df = pd.read_csv(args.input_file, index_col="structure")
    X = df[['PC1', 'PC2']].values

    inertias = []
    silhouettes = []

    Ks = range(2, args.max_clusters + 1)
    for k in Ks:
        kmeans = KMeans(n_clusters=k, random_state=0).fit(X)
        inertias.append(kmeans.inertia_)
        labels = kmeans.labels_
        silhouettes.append(silhouette_score(X, labels))

    # ---- Elbow Plot ----
    plt.figure(figsize=(6, 5))
    plt.plot(Ks, inertias, 'o-', color='blue')
    plt.xlabel('Number of Clusters (k)')
    plt.ylabel('Inertia (Within-cluster Sum of Squares)')
    plt.title(f'Elbow Method for {args.tag}')
    plt.grid(True)
    plt.savefig(f"{args.tag}-elbow.png")

    # ---- Silhouette Plot ----
    plt.figure(figsize=(6, 5))
    plt.plot(Ks, silhouettes, 'o-', color='green')
    plt.xlabel('Number of Clusters (k)')
    plt.ylabel('Average Silhouette Score')
    plt.title(f'Silhouette Method for {args.tag}')
    plt.grid(True)
    plt.savefig(f"{args.tag}-silhouette.png")

    # ---- Print summary ----
    print("\nK  Inertia  Silhouette")
    print("-----------------------")
    for k, inertia, sil in zip(Ks, inertias, silhouettes):
        print(f"{k:<3d} {inertia:<10.2f} {sil:.4f}")

    best_k = Ks[silhouettes.index(max(silhouettes))]
    print(f"\nBest k based on silhouette score = {best_k}")
    print(f"Plots saved as {args.tag}-elbow.[png] and {args.tag}-silhouette.[png]")

    print(f"{best_k}", flush=True)

