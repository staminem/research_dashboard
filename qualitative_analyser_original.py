# -*- coding: utf-8 -*-
"""
NAME:            qualitative_analyser.py
AUTHOR:          Alan Davies
EMAIL:           alan.davies-2@manchester.ac.uk
PROFILE:         https://www.research.manchester.ac.uk/portal/alan.davies-2.html | https://alandavies.netlify.com/
DATE:            12/08/2020
INSTITUTION:     University of Manchester, School of Health Sciences
DESCRIPTION:     Prototype script for formative assessment and feedback of qualitative assigments. 
"""
from tkinter import filedialog as fd
from QualitativeAnalyser import FileLoader, CourseAssignment, AssignmentAnalysis, ScoreAssignments

# Ask user to select the xml file containing assessment data
filename = fd.askopenfilename()

# Create an assignment object and output the details
MIE_course_assignment = CourseAssignment(filename)
MIE_course_assignment.output_assessment_details()

# Ask user to select a folder containing the student assignments
student_assignment_foldername = fd.askdirectory()

# Load the students assignments
reflection_files_loader = FileLoader(student_assignment_foldername)
reflection_files = reflection_files_loader.get_files()

# Add the individual student assignment to the course assingment object
MIE_course_assignment.add_student_assignments(reflection_files)

# Return the students assigments 
assignments = MIE_course_assignment.get_student_assignments()
keywords = MIE_course_assignment.get_keywords()
applies_to = MIE_course_assignment.get_applies_to()

# Create an analyzer for the assigment and run analysis
assignement_analyzer = AssignmentAnalysis()

# Loop over all the assignments
for assignment in assignments:
    assignemnt_text = ' '.join(assignment.paragraph_list)
    assignemnt_text = assignemnt_text.lower()
    matched_contex = assignement_analyzer.keyword_context(assignemnt_text, applies_to)
    intext_refs = assignment.identify_references_in_text(assignemnt_text)
    word_count = assignment.get_word_count(assignemnt_text)
    tokenized_text = assignement_analyzer.tokenize_text(assignemnt_text)
    tokenized_text_no_sw = assignement_analyzer.remove_stop_words(tokenized_text)
    text_as_string = ' '.join(tokenized_text_no_sw)
    nss = assignement_analyzer.critical_analysis(text_as_string)
    main_topics = assignement_analyzer.principle_topics([tokenized_text_no_sw])
    removed_duplicates = assignement_analyzer.remove_duplicate_words(tokenized_text_no_sw)
    pc_keywords = assignement_analyzer.keyword_analysis(keywords, removed_duplicates, keyword_ratio=80)

    # Store the assignments raw scores
    results = {'keywords': pc_keywords,
               'matched_context': matched_contex,
               'topics': main_topics,
               'critical_analysis': nss,
               'wordcount': word_count}

    score_assignments = ScoreAssignments(assignment, MIE_course_assignment)
    score_assignments.compute_assignment_score(results)
    score_assignments.scores_results_table()
    score_assignments.show_feedback(results)
