import pytest 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from main import filter_fastq

def test_filter_by_gc(tmp_path):
    records = [SeqRecord(Seq("ACGT" * 10), id="high_gc", letter_annotations={"phred_quality": [30] * 40}),
                        SeqRecord(Seq("ATAT" * 10), id="low_gc", letter_annotations={"phred_quality": [30] * 40}),
                        SeqRecord(Seq("AGCT" * 10), id="medium_gc", letter_annotations={"phred_quality": [30] * 40})]
    
    input_file = tmp_path / "input.fastq"
    SeqIO.write(records, input_file, "fastq")
    
    output_file = tmp_path / "output.fastq"
    
    filter_fastq(
        input_fastq=str(input_file),
        output_fastq=str(output_file),
        gc_bounds=(50, 100),        
        length_bounds=(0, 100),
        quality_threshold=0)
    
    result = list(SeqIO.parse(output_file, "fastq"))
    
    assert len(result) == 2
    assert result[0].id == "high_gc"
    
def test_filter_by_length(tmp_path):
    records = [SeqRecord(Seq("ACGT" * 10), id="long", letter_annotations={"phred_quality": [30] * 40}),
                        SeqRecord(Seq("ATAT"), id="short", letter_annotations={"phred_quality": [30] * 4}),
                        SeqRecord(Seq("AGCT" * 5), id="medium", letter_annotations={"phred_quality": [30] * 20})]
    
    input_file = tmp_path / "input.fastq"
    SeqIO.write(records, input_file, "fastq")
    
    output_file = tmp_path / "output.fastq"
    
    filter_fastq(
        input_fastq=str(input_file),
        output_fastq=str(output_file),
        gc_bounds=(0, 100),
        length_bounds=(10, 30),     
        quality_threshold=0)
    
    result = list(SeqIO.parse(output_file, "fastq"))
    
    assert len(result) == 1
    assert result[0].id == "medium"
    
def test_quality_filter(tmp_path):
    records = [SeqRecord(Seq("ACGT" * 10), id="high_quality", letter_annotations={"phred_quality": [30] * 40}),
                        SeqRecord(Seq("ATAT" * 10), id="low_quality", letter_annotations={"phred_quality": [10] * 40}),
                        SeqRecord(Seq("AGCT" * 10), id="medium_quality", letter_annotations={"phred_quality": [20] * 40})]
    
    input_file = tmp_path / "input.fastq"
    SeqIO.write(records, input_file, "fastq")
    
    output_file = tmp_path / "output.fastq"
    
    filter_fastq(
        input_fastq=str(input_file),
        output_fastq=str(output_file),
        gc_bounds=(0, 100),
        length_bounds=(0, 100),
        quality_threshold=25)   
    
    result = list(SeqIO.parse(output_file, "fastq"))
    
    assert len(result) == 1
    assert result[0].id == "high_quality"
    
def test_empty_input_file(tmp_path):
    input_file = tmp_path / "empty.fastq"
    input_file.touch() 
    
    output_file = tmp_path / "output.fastq"
    
    filter_fastq(
        input_fastq=str(input_file),
        output_fastq=str(output_file),
        gc_bounds=(0, 100),
        length_bounds=(0, 100),
        quality_threshold=0)
    
    result = list(SeqIO.parse(output_file, "fastq"))
    
    assert len(result) == 0
    
def test_mkdir(tmp_path):
    records = [SeqRecord(Seq("ACGT" * 10), id="test", letter_annotations={"phred_quality": [30] * 40})]
    
    input_file = tmp_path / "input.fastq"
    SeqIO.write(records, input_file, "fastq")
    
    output_dir = tmp_path / "output_dir"
    output_file = output_dir / "output.fastq"
    
    filter_fastq(
        input_fastq=str(input_file),
        output_fastq=str(output_file),
        gc_bounds=(0, 100),
        length_bounds=(0, 100),
        quality_threshold=0)
    
    assert output_file.exists()

def test_nonexistent_input_file(tmp_path):
    output_file = tmp_path / "output.fastq"
    
    with pytest.raises(FileNotFoundError):
        filter_fastq(
            input_fastq=str(tmp_path / "nonexistent.fastq"),
            output_fastq=str(output_file),
            gc_bounds=(0, 100),
            length_bounds=(0, 100),
            quality_threshold=0)
        
def test_combination(tmp_path):
    records = [SeqRecord(Seq("ACGT" * 5), id="high_gc_long_high_quality", letter_annotations={"phred_quality": [30] * 20}),
                        SeqRecord(Seq("ATAT"), id="low_gc_short_low_quality", letter_annotations={"phred_quality": [10] * 4}),
                        SeqRecord(Seq("AGCT" * 5), id="medium_gc_medium_length_medium_quality", letter_annotations={"phred_quality": [20] * 20})]
    
    input_file = tmp_path / "input.fastq"
    SeqIO.write(records, input_file, "fastq")
    
    output_file = tmp_path / "output.fastq"
    
    filter_fastq(
        input_fastq=str(input_file),
        output_fastq=str(output_file),
        gc_bounds=(50, 100),        
        length_bounds=(10, 30),     
        quality_threshold=25)   
    
    result = list(SeqIO.parse(output_file, "fastq"))
    
    assert len(result) == 1
    assert result[0].id == "high_gc_long_high_quality"
    
def test_invalid_gc_bounds(tmp_path):
    records = [SeqRecord(Seq("ACGT" * 10), id="test", letter_annotations={"phred_quality": [30] * 40})]
    
    input_file = tmp_path / "input.fastq"
    SeqIO.write(records, input_file, "fastq")
    
    output_file = tmp_path / "output.fastq"
    
    with pytest.raises(ValueError):
        filter_fastq(
            input_fastq=str(input_file),
            output_fastq=str(output_file),
            gc_bounds=(100, 50),        
            length_bounds=(0, 100),     
            quality_threshold=0)