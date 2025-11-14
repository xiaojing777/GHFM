# GHFM: Generalized Heterogeneous Functional Model

This repository contains R implementations for the heterogeneous functional regression framework described in the manuscript. The code follows the workflow below.

## Workflow

1. **Start (GHFM)**
2. **Check sample size**
   - **Large sample**: derive pre-clustering groups.
   - **Not large**: skip pre-clustering.
3. **Check response type**
   - **Continuous**: use HFGM.
   - **Binary**: use HFLM.
4. **Fit model**
   - If pre-clustering was used:
     - Continuous → apply HFGM to the \(K\) pre-clustering groups.
     - Binary → apply HFLM to the \(K\) pre-clustering groups.
   - If pre-clustering was not used:
     - Continuous → apply HFGM to all subjects.
     - Binary → apply HFLM to all subjects.
5. **(Optional) Testing**
   - Use the testing script to compare the group-varying model vs. the pooled model.

The flow above corresponds to Figure 2 in the manuscript.

## Files

- `pre_clustering.R`  
  - Functions to obtain \(K\) initial subject groups for large-\(n\) cases.
  - Supports Gaussian and logistic responses.

- `HFGM_LQA.R`  
  - LQA-based implementation of the heterogeneous **functional** regression for continuous response.
  - Use this when there is **no** pre-clustering or when you want a direct LQA version.

- `HFLM_LQA.R`  
  - LQA-based implementation of the heterogeneous **functional logistic** regression for binary response.

- `HFGM_ADMM.R`  
  - ADMM-based implementation for continuous response.
  - Contains two versions: without pre-groups and with pre-groups.

- `HFLM_ADMM.R`  
  - ADMM-based implementation for binary response with pre-clustering groups.

- `GHFM_testing.R`  
  - Procedures for comparing the reduced model (common functional coefficient) vs. the full model (group-specific functional coefficients) for both Gaussian and Bernoulli responses.

## Typical usage

- **Large sample, continuous response**
  1. Run `pre_clustering.R` to get \(K\) groups.
  2. Run either `HFGM_LQA.R` or `HFGM_ADMM.R` on these \(K\) groups.

- **Large sample, binary response**
  1. Run `pre_clustering.R` to get \(K\) groups.
  2. Run `HFLM_LQA.R` or `HFLM_ADMM.R` on these \(K\) groups.

- **Small/medium sample**
  - Skip pre-clustering, run the corresponding LQA/ADMM file on all subjects.

- **Model comparison**
  - After fitting, use `GHFM_testing.R` to test group-varying model vs. pooled model.

## Notes

- This README assumes the tuning parameters (smoothing, fusion, ADMM penalties) have been selected.
- All scripts are written in base R style and use the `fda` package for functional data objects.

## Simulation example

A simple end-to-end example is provided in `Simulation_Example.R`.

After cloning this repository, you can run it from the GHFM root directory via:

```bash
Rscript Simulation_Example.R


