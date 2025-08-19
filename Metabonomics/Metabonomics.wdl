# Major steps:
#  1) (Optional) RAW -> mzML (ProteoWizard msconvert)
#  2) Peak picking (OpenMS PeakPickerHiRes) if profile data
#  3) Feature detection (OpenMS FeatureFinderMetabo)
#  4) Retention-time alignment (OpenMS MapAlignerPoseClustering)
#  5) Adduct deconvolution (OpenMS MetaboliteAdductDecharger)
#  6) Feature linking across samples (OpenMS FeatureLinkerUnlabeledKD)
#  7) Export to feature matrix (OpenMS TextExporter)
#  8) (Optional) Blank subtraction and simple QC summary
#
# Notes:
#  - Uses Docker images: chambm/pwiz-skyline-i-agree-to-the-vendor-licenses (msconvert)
#                        openms/openms:latest
#  - Assumes centroided data downstream of step 2. Set do_peak_picking=false if input mzMLs are already centroided.
#  - Provide either vendor RAWs (e.g., Thermo .raw) or mzMLs. If both are provided, mzMLs take precedence and RAW->mzML is skipped.
################################################################################

workflow MetabolomicsPipeline {
  meta {
    description: "End-to-end upstream metabolomics LC–MS pipeline using OpenMS tools."
  }

  input {
    # One of the following two should be provided
    Array[File]? vendor_raws                               # e.g., .raw files from Thermo
    Array[File]? mzml_inputs                               # Already converted mzMLs

    # Sample sheet (2 columns: sample_id, file_basename) is optional but recommended
    File? sample_sheet_tsv

    # Processing toggles
    Boolean do_peak_picking = true                         # Apply PeakPickerHiRes
    Boolean enable_blank_subtraction = false               # Requires blanks list

    # Optional blanks/QCs
    Array[String]? blank_sample_ids
    Array[String]? qc_sample_ids

    # PeakPickerHiRes params
    Float pphr_signal_to_noise = 0.0                       # 0 = OpenMS default; set >3 for noisy data

    # FeatureFinderMetabo params (tune as appropriate)
    Float ffm_noise_threshold_int = 1000.0
    Float ffm_mass_error_ppm = 10.0
    Float ffm_min_trace_length = 2.0

    # MapAligner params
    Float mapalign_max_rt_shift = 120.0                    # seconds

    # FeatureLinker params
    Float fl_tolerance_mz = 0.005
    Float fl_tolerance_rt = 5.0                            # seconds

    # Resources
    Int cpus = 4
    String ram_gb = "8G"
    String disk_gb = "50G"
  }

  ############################
  # Select inputs (RAW->mzML)
  ############################
  if (defined(mzml_inputs)) {
    call IdentityMzML as useProvidedMzML { input: mzmls = select_first([mzml_inputs, []]) }
  } else {
    if (defined(vendor_raws)) {
      scatter (raw in select_first([vendor_raws, []])) {
        call Msconvert as convert { input: raw_file = raw, cpus = cpus, disk_gb = disk_gb }
      }
      call GatherMzML { input: mzmls = convert.mzml }
      call IdentityMzML as useProvidedMzML { input: mzmls = GatherMzML.mzmls }
    } else {
      call FailNoInputs
    }
  }

  ############################################
  # Optional Peak Picking (profile -> centroid)
  ############################################
  Array[File] centroided_mzmls
  if (do_peak_picking) {
    scatter (f in useProvidedMzML.mzmls) {
      call PeakPickerHiResTask as pphr {
        input:
          mzml = f,
          signal_to_noise = pphr_signal_to_noise,
          cpus = cpus,
          ram_gb = ram_gb,
          disk_gb = disk_gb
      }
    }
    call GatherMzML as gatherCentroided { input: mzmls = pphr.centroided_mzml }
    centroided_mzmls = gatherCentroided.mzmls
  } else {
    centroided_mzmls = useProvidedMzML.mzmls
  }

  ######################
  # Feature Detection
  ######################
  scatter (f in centroided_mzmls) {
    call FeatureFinderMetaboTask as ffm {
      input:
        mzml = f,
        noise_threshold_int = ffm_noise_threshold_int,
        mass_error_ppm = ffm_mass_error_ppm,
        min_trace_length = ffm_min_trace_length,
        cpus = cpus,
        ram_gb = ram_gb,
        disk_gb = disk_gb
    }
  }
  call GatherFeatureXML { input: features = ffm.featurexml }

  ######################
  # RT Alignment
  ######################
  call MapAlignerPoseClusteringTask as align {
    input:
      featurexml_list = GatherFeatureXML.features,
      max_rt_shift = mapalign_max_rt_shift,
      cpus = cpus,
      ram_gb = ram_gb,
      disk_gb = disk_gb
  }

  #################################
  # Adduct Deconvolution (optional)
  #################################
  scatter (fx in align.aligned_featurexml) {
    call MetaboliteAdductDechargerTask as decharge {
      input:
        featurexml = fx,
        cpus = cpus,
        ram_gb = ram_gb,
        disk_gb = disk_gb
    }
  }
  call GatherFeatureXML as gatherDecharged { input: features = decharge.decharged_featurexml }

  ###############################
  # Feature Linking across samples
  ###############################
  call FeatureLinkerUnlabeledKDTask as link {
    input:
      featurexml_list = gatherDecharged.features,
      mz_tolerance = fl_tolerance_mz,
      rt_tolerance = fl_tolerance_rt,
      cpus = cpus,
      ram_gb = ram_gb,
      disk_gb = disk_gb
  }

  ########################
  # Export feature matrix
  ########################
  call TextExporterTask as export {
    input:
      consensusxml = link.consensusxml,
      cpus = cpus,
      ram_gb = ram_gb,
      disk_gb = disk_gb
  }

  ############################
  # Optional: Blank Subtraction
  ############################
  File? blanks_removed_tsv
  if (enable_blank_subtraction && defined(sample_sheet_tsv) && defined(blank_sample_ids)) {
    call BlankSubtractionTask as blanks {
      input:
        feature_tsv = export.feature_tsv,
        sample_sheet_tsv = select_first([sample_sheet_tsv, ""]),
        blank_ids = select_first([blank_sample_ids, []])
    }
    blanks_removed_tsv = blanks.filtered_feature_tsv
  }

  output {
    Array[File] final_mzmls = centroided_mzmls
    Array[File] featurexml_per_sample = gatherDecharged.features
    File aligned_consensusxml = link.consensusxml
    File feature_matrix_tsv = export.feature_tsv
    File? feature_matrix_blanks_removed_tsv = blanks_removed_tsv
    File? map_alignment_tsv = align.trafo_info_tsv
    File? qc_summary_tsv = blanks.qc_summary_tsv
  }
}

################################################################################
# Tasks
################################################################################

task Msconvert {
  input {
    File raw_file
    Int cpus = 4
    String disk_gb = "30G"
  }
  command <<<
    set -euo pipefail
    mkdir -p out
    base=$(basename ~{raw_file})
    stem=${base%.*}
    msconvert "~{raw_file}" \
      --mzML \
      --zlib \
      --filter "peakPicking true 1-" \
      --outfile "${stem}.mzML" \
      --outdir out
    echo "Converted: out/${stem}.mzML"
  >>>
  output {
    File mzml = "out/*.mzML"
  }
  runtime {
    docker: "chambm/pwiz-skyline-i-agree-to-the-vendor-licenses"
    cpu: cpus
    disks: "local-disk ~{disk_gb}"
  }
}

# Pass-through to unify interface
task IdentityMzML {
  input { Array[File] mzmls }
  command <<<
    ls -1 ~{sep=' ' mzmls} > mzml.list
  >>>
  output {
    Array[File] mzmls = mzmls
  }
  runtime { docker: "openms/openms:latest" }
}

task GatherMzML {
  input { Array[File] mzmls }
  command <<<
    ls -1 ~{sep=' ' mzmls} > gathered_mzmls.list
  >>>
  output { Array[File] mzmls = mzmls }
  runtime { docker: "openms/openms:latest" }
}


task PeakPickerHiResTask {
  input {
    File mzml
    Float signal_to_noise = 0.0
    Int cpus = 4
    String ram_gb = "8G"
    String disk_gb = "20G"
  }
  command <<<
    set -euo pipefail
    base=$(basename ~{mzml})
    stem=${base%.mzML}
    PeakPickerHiRes -in "~{mzml}" -out "${stem}.centroid.mzML" -signal_to_noise ~{signal_to_noise}
  >>>
  output { File centroided_mzml = "*.centroid.mzML" }
  runtime {
    docker: "openms/openms:latest"
    cpu: cpus
    memory: ram_gb
    disks: "local-disk ~{disk_gb}"
  }
}


task FeatureFinderMetaboTask {
  input {
    File mzml
    Float noise_threshold_int = 1000.0
    Float mass_error_ppm = 10.0
    Float min_trace_length = 2.0
    Int cpus = 4
    String ram_gb = "8G"
    String disk_gb = "20G"
  }
  command <<<
    set -euo pipefail
    base=$(basename ~{mzml})
    stem=${base%.mzML}
    FeatureFinderMetabo \
      -in "~{mzml}" \
      -out "${stem}.featureXML" \
      -algorithm:common:noise_threshold_int ~{noise_threshold_int} \
      -algorithm:common:mass_error_ppm ~{mass_error_ppm} \
      -algorithm:detect:chrom_peak_snr 3 \
      -algorithm:detect:min_trace_length ~{min_trace_length}
  >>>
  output { File featurexml = "*.featureXML" }
  runtime {
    docker: "openms/openms:latest"
    cpu: cpus
    memory: ram_gb
    disks: "local-disk ~{disk_gb}"
  }
}


task GatherFeatureXML {
  input { Array[File] features }
  command <<<
    ls -1 ~{sep=' ' features} > featurexml.list
  >>>
  output { Array[File] features = features }
  runtime { docker: "openms/openms:latest" }
}


task MapAlignerPoseClusteringTask {
  input {
    Array[File] featurexml_list
    Float max_rt_shift = 120.0
    Int cpus = 4
    String ram_gb = "8G"
    String disk_gb = "20G"
  }
  command <<<
    set -euo pipefail
    mkdir -p aligned
    inputs=(~{sep=' ' featurexml_list})
    # Use the first file as reference
    ref="${inputs[0]}"
    for fx in "${inputs[@]}"; do
      base=$(basename "$fx")
      MapAlignerPoseClustering \
        -in "$fx" \
        -out "aligned/${base%.featureXML}.aligned.featureXML" \
        -reference:file "$ref" \
        -max_num_peaks_considered 3000 \
        -max_rt_shift ~{max_rt_shift}
    done
    # Optionally export transformation info
    MapRTTransformer -in "$ref" -out "aligned/ref.trafoXML" || true
    TrafoXML2TSV -in "aligned/ref.trafoXML" -out "aligned/trafo.tsv" || true
  >>>
  output {
    Array[File] aligned_featurexml = glob("aligned/*.aligned.featureXML")
    File? trafo_info_tsv = "aligned/trafo.tsv"
  }
  runtime {
    docker: "openms/openms:latest"
    cpu: cpus
    memory: ram_gb
    disks: "local-disk ~{disk_gb}"
  }
}


task MetaboliteAdductDechargerTask {
  input {
    File featurexml
    Int cpus = 2
    String ram_gb = "4G"
    String disk_gb = "10G"
  }
  command <<<
    set -euo pipefail
    base=$(basename ~{featurexml})
    stem=${base%.featureXML}
    MetaboliteAdductDecharger -in "~{featurexml}" -out "${stem}.decharged.featureXML" -potential_adducts "[M+H]+,[M+Na]+,[M+K]+,[M-H]-,[M+Cl]-" -negative_em_removal true
  >>>
  output { File decharged_featurexml = "*.decharged.featureXML" }
  runtime {
    docker: "openms/openms:latest"
    cpu: cpus
    memory: ram_gb
    disks: "local-disk ~{disk_gb}"
  }
}


task FeatureLinkerUnlabeledKDTask {
  input {
    Array[File] featurexml_list
    Float mz_tolerance = 0.005
    Float rt_tolerance = 5.0
    Int cpus = 4
    String ram_gb = "8G"
    String disk_gb = "20G"
  }
  command <<<
    set -euo pipefail
    FeatureLinkerUnlabeledKD \
      -in ~{sep=',' featurexml_list} \
      -out linked.consensusXML \
      -algorithm:distance_MZ:unit mz \
      -algorithm:distance_MZ:threshold ~{mz_tolerance} \
      -algorithm:distance_RT:unit sec \
      -algorithm:distance_RT:threshold ~{rt_tolerance}
  >>>
  output { File consensusxml = "linked.consensusXML" }
  runtime {
    docker: "openms/openms:latest"
    cpu: cpus
    memory: ram_gb
    disks: "local-disk ~{disk_gb}"
  }
}


task TextExporterTask {
  input {
    File consensusxml
    Int cpus = 1
    String ram_gb = "2G"
    String disk_gb = "5G"
  }
  command <<<
    set -euo pipefail
    TextExporter -in "~{consensusxml}" -out "features.tsv" -export:consensus:matrix true -export:consensus:decharged_masses true -export:consensus:rt true -export:consensus:charge true -export:consensus:ids true
  >>>
  output { File feature_tsv = "features.tsv" }
  runtime {
    docker: "openms/openms:latest"
    cpu: cpus
    memory: ram_gb
    disks: "local-disk ~{disk_gb}"
  }
}


task BlankSubtractionTask {
  input {
    File feature_tsv
    File sample_sheet_tsv
    Array[String] blank_ids
  }
  command <<<
    set -euo pipefail
    python3 << 'PY'
import csv, sys
import statistics as stats
from collections import defaultdict

feature_file = "~{feature_tsv}"
sample_sheet = "~{sample_sheet_tsv}"
blank_ids = set(~{sep=',' blank_ids})

# Read sample order (map columns)
with open(sample_sheet, 'r', newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    sample_ids = [r['sample_id'] for r in reader]

# Load feature table
with open(feature_file, 'r', newline='') as f:
    rows = list(csv.reader(f, delimiter='\t'))
header = rows[0]
# Identify intensity columns by sample ids
col_idx = {sid: header.index(sid) for sid in sample_ids if sid in header}
blank_cols = [col_idx[s] for s in sample_ids if s in blank_ids and s in col_idx]

filtered = [header]
for r in rows[1:]:
    try:
        blank_vals = [float(r[i]) for i in blank_cols]
        blank_med = stats.median(blank_vals) if blank_vals else 0.0
        # Keep if any non-blank sample exceeds 3x blank median
        keep = False
        for sid, idx in col_idx.items():
            if sid in blank_ids: 
                continue
            try:
                val = float(r[idx])
            except:
                val = 0.0
            if blank_med == 0 and val > 0:
                keep = True; break
            if blank_med > 0 and val > 3 * blank_med:
                keep = True; break
        if keep:
            filtered.append(r)
    except Exception:
        # keep row if parsing fails
        filtered.append(r)

with open('features.blanks_removed.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerows(filtered)

# Simple QC summary
qc = []
for sid, idx in col_idx.items():
    if sid in blank_ids:
        continue
    nonzeros = 0
    zeros = 0
    for r in filtered[1:]:
        try:
            v = float(r[idx])
            if v > 0: nonzeros += 1
            else: zeros += 1
        except:
            zeros += 1
    qc.append((sid, nonzeros, zeros))

with open('qc_summary.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['sample_id','nonzero_features','zero_features'])
    for sid, n, z in qc:
        w.writerow([sid, n, z])
PY
  >>>
  output {
    File filtered_feature_tsv = "features.blanks_removed.tsv"
    File qc_summary_tsv = "qc_summary.tsv"
  }
  runtime { docker: "python:3.11-slim" }
}


task FailNoInputs {
  command <<<
    echo "ERROR: Provide either 'vendor_raws' or 'mzml_inputs'" >&2
    exit 1
  >>>
  output { String msg = read_string(stdout()) }
  runtime { docker: "alpine:3.20" }
}

################################################################################
# Parameter Metadata (helps when UIs render forms)
################################################################################
parameter_meta {
  MetabolomicsPipeline.vendor_raws: "Vendor RAW files (e.g., Thermo .raw). Use when mzMLs are not provided."
  MetabolomicsPipeline.mzml_inputs: "Pre-converted centroided mzML files. Takes precedence over vendor_raws."
  MetabolomicsPipeline.sample_sheet_tsv: "TSV with columns: sample_id, file_basename (basename without extension)."
  MetabolomicsPipeline.do_peak_picking: "Run PeakPickerHiRes if data are profile-mode."
  MetabolomicsPipeline.enable_blank_subtraction: "Enable simple blank-based feature filtering."
  MetabolomicsPipeline.blank_sample_ids: "Sample IDs that are blanks (must match header columns in features.tsv)."
  MetabolomicsPipeline.qc_sample_ids: "Optional pooled QC sample IDs (not used directly; reserved for future QC steps)."
  MetabolomicsPipeline.pphr_signal_to_noise: "PeakPickerHiRes S/N threshold (0 = default)."
  MetabolomicsPipeline.ffm_noise_threshold_int: "FeatureFinderMetabo intensity noise threshold."
  MetabolomicsPipeline.ffm_mass_error_ppm: "FeatureFinderMetabo mass error in ppm."
  MetabolomicsPipeline.ffm_min_trace_length: "FeatureFinderMetabo minimum chromatographic trace length."
  MetabolomicsPipeline.mapalign_max_rt_shift: "Maximum allowed RT shift (sec) during alignment."
  MetabolomicsPipeline.fl_tolerance_mz: "Feature linking m/z tolerance (in Da)."
  MetabolomicsPipeline.fl_tolerance_rt: "Feature linking RT tolerance (sec)."
}

################################################################################
# Example inputs JSON (save as inputs.json)
################################################################################
# {
#   "MetabolomicsPipeline.mzml_inputs": [
#     "gs://your-bucket/run1/sample01.mzML",
#     "gs://your-bucket/run1/sample02.mzML",
#     "gs://your-bucket/run1/blank01.mzML"
#   ],
#   "MetabolomicsPipeline.sample_sheet_tsv": "gs://your-bucket/meta/samples.tsv",
#   "MetabolomicsPipeline.do_peak_picking": false,
#   "MetabolomicsPipeline.enable_blank_subtraction": true,
#   "MetabolomicsPipeline.blank_sample_ids": ["blank01"],
#   "MetabolomicsPipeline.ffm_noise_threshold_int": 1000.0,
#   "MetabolomicsPipeline.ffm_mass_error_ppm": 10.0,
#   "MetabolomicsPipeline.ffm_min_trace_length": 2.0,
#   "MetabolomicsPipeline.fl_tolerance_mz": 0.005,
#   "MetabolomicsPipeline.fl_tolerance_rt": 5.0,
#   "MetabolomicsPipeline.cpus": 8,
#   "MetabolomicsPipeline.ram_gb": "16G",
#   "MetabolomicsPipeline.disk_gb": "100G"
# }

################################################################################
# Usage (Cromwell example):
#   java -jar cromwell.jar run metabolomics.wdl -i inputs.json
################################################################################
