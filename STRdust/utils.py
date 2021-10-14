import os
import sys
import logging
import tempfile
import pysam


class InputException(Exception):
    pass


def setup_temp_dirs(out_dir):
    ins_dir = os.path.join(out_dir, "chrs_ins_tmp")
    if not os.path.isdir(ins_dir):
        os.mkdir(ins_dir)

    vcf_dir = os.path.join(out_dir, "chrs_vcf_tmp")
    if not os.path.isdir(vcf_dir):
        os.mkdir(vcf_dir)
    return ins_dir, vcf_dir


def create_output_directory(out_dir):
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    else:
        path_to_dir = os.path.join(out_dir, "test")
        if not _validate_path(path_to_dir):
            sys.exit(f"Problem with writing permissions in output directory. {path_to_dir}\n")
    return os.path.abspath(out_dir)


def _validate_path(path):
    try:
        os.makedirs(path, exist_ok=True)
        temp_dir_path = tempfile.mkdtemp(dir=path)
        os.rmdir(temp_dir_path)
        return True
    except OSError:
        return False


def _enable_logging(out_dir, debug, overwrite):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    logger = logging.getLogger()
    log_file = os.path.join(out_dir, "STRdust.log")
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%Y-%m-%d %H:%M:%S")

    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)

    if not debug:
        console_log.setLevel(logging.INFO)

    if overwrite:
        open(log_file, "w").close()

    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)


def _check_bam_files(bam_file):
    """
    Check existance of input files and generate index file if it is absent
    :param bam_file: phased bam file with/without bai file
    """

    if not os.path.exists(bam_file):
        raise InputException(f"Can't open {bam_file}")

    samfile = pysam.AlignmentFile(bam_file, "rb")
    if not samfile.has_index():
        logging.info("Input bam file does not have index file (.bai). Generating now.")
        pysam.index(bam_file)
