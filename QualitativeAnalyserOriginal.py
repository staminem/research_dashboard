"""
NAME:            QualitativeAnalyser.py
AUTHOR:          Alan Davies
EMAIL:           alan.davies-2@manchester.ac.uk
PROFILE:         https://www.research.manchester.ac.uk/portal/alan.davies-2.html | https://alandavies.netlify.com/
DATE:            12/08/2020
INSTITUTION:     University of Manchester, School of Health Sciences
DESCRIPTION:     Functionality of importing and analysing student assignment.
"""

import ntpath
import gensim
import string
import nltk
import os, sys, re
from gensim import corpora
from docx import Document
from fuzzywuzzy import fuzz
from nltk.tokenize import word_tokenize 
from nltk.corpus import stopwords
from nltk.sentiment.vader import SentimentIntensityAnalyzer
import xml.etree.ElementTree as et
from tabulate import tabulate 

class FileLoader:
    """Class to manage importing of text files in *.docx format."""    
    def __init__(self, folder_selected):
        """Initialization function."""
        self.folder_selected = folder_selected
        self.file_list = []
        
    def get_files(self):
        """Returns the files as a list."""
        for path, subdirs, files in os.walk(self.folder_selected):
            for name in files:
                if name.endswith('.docx'):
                    self.file_list.append(os.path.join(path, name))                    
        return self.file_list

class ParagraphParser:
    """Add document paragraphs to list object."""
    def __init__(self):
        """Initialization function."""
        self.document_paragraphs = []
        
    def get_paragraphs(self, document_list):
        """Return the documents and associated paragraphs."""
        for file in document_list:
            doc = Document(file)
            for paragraph in doc.paragraphs:
                self.document_paragraphs.append(paragraph.text)                                
        return self.document_paragraphs

class CourseAssignment:
    """Represents a course unit assignement."""
    def __init__(self, assigment_data_file_xml):
        """Initialization function."""
        self.student_assignments = []
        self.marks_per_concept = []
        self.keywords = []
        self.applies_to = []
        self.applies_to_concept = []
        self.applies_list = []
        self.feedback = []
        self.filename = assigment_data_file_xml
        try:
            self.content = et.parse(self.filename)
            self.root = self.content.getroot()
        except:
            print("Error parsing file:", self.filename)
            sys.exit(1)
             
        self.course_name = self.root.find('course')
        self.course_code = self.root.find('code')
        self.assessment_type = self.root.find('type')
        self.assessment_description = self.root.find('description')
        self.number_of_concepts = len(self.root.findall('concept'))            
        self.topics = self.root.find('topics')        
        self.word_count = self.root.find('wordcount')
        self.word_count_limit = int(self.word_count.text)
        self.word_count_penalty_type = self.word_count.attrib['penalty']
        self.word_count_include_refs = self.word_count.attrib['include_refs']
        self.actual_word_count = 0
        self.concept_titles = []
        v = []  
        #concept_count = 1
        concept_feedback = []
        for concept in self.root.iter('concept'):                        
            self.marks_per_concept.append(float(concept.attrib['max_marks']))            
            self.keywords.append([i.text for i in concept.iter('keyword')])      
            self.concept_titles.append([i.text for i in concept.iter('title')])
            for fb in concept.iter('feedback'):
                vals = fb.attrib.values()
                vals_list = list(vals)
                if '-' in vals_list[0]:
                    x = vals_list[0].split('-')
                    mark = list(range(int(x[0]), int(x[1])+1))
                else:
                    mark = [int(vals_list[0])]
                concept_feedback.append({'mark': mark, 'feedback_text': fb.text})
            self.feedback.append(concept_feedback) 
            concept_feedback = []
            #concept_count+=1
            for ap in concept.iter('applies'):
                vals = ap.attrib.values()                
                vals_list = list(vals)
                v.append(vals_list[0])
                st_list = [i.lstrip() for i in ap.text.split(',')]                
                self.applies_to.append({vals_list[0]: tuple(st_list)})     
            self.applies_list.append(v)  
            v = []  
            self.applies_to_concept.append(self.applies_to) 
            self.applies_to = []                            
        
    def output_assessment_details(self):
        """Output meta data about the assessment."""
        print("Assessment details:")
        print("Course name:", self.course_name.text)
        print("Course code:", self.course_code.text)
        print("Assessment type:", self.assessment_type.text)
        print("Assessment description:", self.assessment_description.text)
        
    def add_student_assignments(self, assignments):
        """Create a list of individual student assignment objects."""
        for assigment_file in assignments:
            current_assignment = StudentAssignment(assigment_file)
            self.student_assignments.append(current_assignment)
        
    def get_student_assignments(self):
        """Return the list of student assignments."""
        return self.student_assignments
    
    def get_keywords(self):
        """Return the list of keywords."""
        return self.keywords
    
    def get_applies_to(self):
        """Return the applies to context."""
        return self.applies_to_concept
    
    def get_applies_keywords(self):
        """Return the list of applies keys."""
        return self.applies_list
        
class StudentAssignment:
    """Represents an individual student assignement."""
    def __init__(self, assigment_file):
        """Initialization function."""
        self.assigment_file = assigment_file
        self.assignment_mark_pc = 0
        self.word_count_penalty = 0
        self.keyword_score = 0
        self.context_score = 0
        self.topic_score = 0
        self.concept_scores = []
        self.critical_analysis_score = 0
        self.student_feedback = ""
        paragraph_parser = ParagraphParser()
        self.paragraph_list = paragraph_parser.get_paragraphs([self.assigment_file])
        self.student_id = self.path_leaf(self.assigment_file)
        
    def path_leaf(self, file_path):
        """Function to return just filename from file path."""
        head, tail = ntpath.split(file_path)
        return tail.split('.',1)[0] or ntpath.basename(head).split('.',1)[0]
    
    def identify_references_in_text(self, assignment_text):
        """Extract the in text references from the assignemnt text."""
        text_refs_search = re.compile(r"\([^()\d]*\d[^()]*\)")
        text_refs = re.findall(text_refs_search, assignment_text)
        text_refs_list = list(text_refs)
        text_refs_count = {i:text_refs_list.count(i) for i in text_refs_list}
        return text_refs_count
    
    def get_word_count(self, assigment_file, minus_references=True):
        """Return the word count with or without references and ref list."""
        remove_characters = 0
        total_words = len(assigment_file.split())
        if not minus_references:
            return total_words
        ref_index = assigment_file.rfind('references')
        if ref_index == -1:
            ref_index = assigment_file.rfind('bibliography')
            if ref_index == -1:
                return -1         
        intext_refs = self.identify_references_in_text(assigment_file)
        for key, value in intext_refs.items():
            remove_characters = (remove_characters + len(key.split())) * value
        remove_characters = remove_characters + len(assigment_file[ref_index: ].split())
        return total_words - remove_characters
    
class AssignmentAnalysis:
    """Class of analysis methods."""
    def __init__(self):
        """Initialization function."""
        self.keyword_fuzz_ratios = []
    
    def tokenize_text(self, assignment_text):
        """Split text into words, remove punctuation."""
        word_tokens = word_tokenize(assignment_text)
        words = [word for word in word_tokens if word.isalpha()]
        return words
    
    def remove_stop_words(self, assignment_text):
        """Remove all stop words form assignment text."""
        stop_words = set(stopwords.words('english'))
        words = [i for i in assignment_text if not i in stop_words]
        return words
    
    def remove_duplicate_words(self, assignment_text):
        """Remove all the duplicate words from text for keyword mathing."""
        filtered_list = []
        for word in assignment_text:
            if word not in filtered_list:
                filtered_list.append(word)
        return filtered_list    
            
    def keyword_analysis(self, keywords, assignment_text, keyword_ratio=80, echo=False):
        """Match keywords using fuzzy strings."""
        word_matched = 0
        concept_matches = []        
        self.keyword_ratio = keyword_ratio
        for keywords_in_concept in keywords:  
            for keyword in keywords_in_concept:                
                for word in assignment_text:  
                    ratio = fuzz.ratio(word, keyword)
                    if ratio >= self.keyword_ratio: 
                        word_matched = word_matched + 1
                        if echo:
                            print("Matched (word/kw): ", word, keyword, ratio)
                        self.keyword_fuzz_ratios.append({'word': word, 'keyword': keyword, 'fuzzy-ratio': ratio})                               
            if word_matched > 0 and word_matched < len(keywords_in_concept): 
                concept_matches.append((word_matched / len(keywords_in_concept)) * 100)
            elif word_matched > 0 and word_matched >= len(keywords_in_concept): 
                concept_matches.append(100)
            elif word_matched == 0:
                concept_matches.append(word_matched)
            word_matched = 0            
        return tuple(concept_matches)
    
    def get_keyword_fuzzy_ratios(self):
        """Return the fuzzy ratio's of the matched keywords."""
        return self.keyword_fuzz_ratios
    
    def keyword_context(self, assignment_text, applies_to, match_ratio=80):
        """Fuzzy match bi-grams."""
        matched_bigrams = []     
        bigrams_per_context = []        
        for applied_context in applies_to:
            assignment_text_no_punctuation = assignment_text.translate(str.maketrans('', '', string.punctuation))         
            nltk_tokens = nltk.word_tokenize(assignment_text_no_punctuation)
            kc_bigrams = list(nltk.bigrams(nltk_tokens))                
            for key_context in applied_context:
                for key, value in key_context.items():
                    for item in value:
                        for gram in kc_bigrams:
                            r1 = fuzz.ratio(key, gram[0])
                            r2 = fuzz.ratio(item, gram[1])
                            if r1 >= match_ratio and r2 >= match_ratio:                            
                                matched_bigrams.append({'term_1': key, 'bigram_1': gram[0], 'term_2': item, 'bigram_2': gram[1], 'ratio_1': r1, 'ratio_2': r2})
            bigrams_per_context.append(matched_bigrams) 
            matched_bigrams = [] 
        return bigrams_per_context
    
    def principle_topics(self, assignment_text, num_topics=5, num_words=4):
        """Uses topic modelling to check main themes of work."""
        dictionary = corpora.Dictionary(assignment_text)
        corpus = [dictionary.doc2bow(text) for text in assignment_text]
        lda_model = gensim.models.ldamodel.LdaModel(corpus, num_topics=num_topics, id2word=dictionary, passes=15)
        topics = lda_model.print_topics(num_words=num_words)
        return topics   
    
    def critical_analysis(self, assignment_text):
        """Runs a sentiment anaylsis with net sentiment as proxy for critical analysis."""
        sent_score = SentimentIntensityAnalyzer().polarity_scores(assignment_text)
        nss = (sent_score['pos'] * 100) - (sent_score['neg'] * 100)
        sent_score.update({'nss': nss})
        return sent_score
        
class ScoreAssignments:
    """Produces scores based on the results of the various analysis methods."""
    def __init__(self, assignment, course_assignment_object):
        """Initialization function."""
        self.student_assignment = assignment
        self.course_assignment_object = course_assignment_object
           
    def compute_assignment_score(self, results, extra_words_pc=10, pc_penalty_per_words=[1, 100]): 
        """Process the results to extract final scores and add to individual assignment."""        
        if 'topics' in results.keys():
            topic_list = []
            topic_terms_pattern = re.compile('"([^"]*)"')
            for topic in results['topics']:                
                topic_list.append(re.findall(topic_terms_pattern, topic[1]))
            flatten_list = [item for sublist in topic_list for item in sublist]                    
            assignment_topic_list = self.course_assignment_object.topics.text.split(',')
            assignment_topic_list = [i.strip(' ') for i in assignment_topic_list]
            count_common_items = set(flatten_list) & set(assignment_topic_list)
            self.student_assignment.topic_score = (len(count_common_items) / len(assignment_topic_list))*100
        if 'critical_analysis' in results.keys():
            nss = results['critical_analysis']['nss']
            nss_scores = list(range(0,105,5))
            s1 = list(range(0,110,10))
            s2 = list(range(90,-10,-10))
            scores = s1 + s2
            rounded_score = 5*round(nss/5)
            idx = nss_scores.index(rounded_score)
            self.critical_analysis_score = scores[idx]            
        if 'matched_context' in results.keys():
            ap_keywords = self.course_assignment_object.get_applies_keywords()
            keys_to_match_per_concept = [len(i) for i in ap_keywords]
            concept_matches = [len(i) for i in results['matched_context']]
            match_pc = []
            if len(keys_to_match_per_concept) == len(concept_matches):
                for i,j in zip(keys_to_match_per_concept, concept_matches):            
                    match_pc.append((j/i)*100)
            else:
                 match_pc.append([0])      
        if 'wordcount' in results.keys():            
            # Add an extra n% of words to the limit
            self.course_assignment_object.actual_word_count = results['wordcount'] 
            word_limit = self.course_assignment_object.word_count_limit + ((self.course_assignment_object.word_count_limit * extra_words_pc) / 100)            
            if results['wordcount'] > word_limit:
                remaining_amount = results['wordcount'] - word_limit                
                self.student_assignment.word_count_penalty = (remaining_amount // pc_penalty_per_words[1])*pc_penalty_per_words[0]
            else:
                self.student_assignment.word_count_penalty = 0                
        for i in range(0, self.course_assignment_object.number_of_concepts):
            if 'keywords' in results.keys():        
                self.student_assignment.concept_scores.append(((results['keywords'][i] + match_pc[i])/200)*100)
        score = (sum(self.student_assignment.concept_scores)/300)*100
        topics_and_analysis_score = (self.critical_analysis_score + self.student_assignment.topic_score) / self.course_assignment_object.number_of_concepts
        self.student_assignment.assignment_mark_pc = (score + topics_and_analysis_score) - self.student_assignment.word_count_penalty
        
    def scores_results_table(self):
        """Output the results in the table.""" 
        headers = [i[0] for i in self.course_assignment_object.concept_titles]
        headers.append('Word count penalty')
        headers.append('Final mark %')
        output_table = [str(round(i)) for i in self.student_assignment.concept_scores]
        output_table.append(str(self.student_assignment.word_count_penalty))
        output_table.append(str(round(self.student_assignment.assignment_mark_pc)))
        print("Student id:",self.student_assignment.student_id, " Word count:",self.course_assignment_object.actual_word_count)
        print(tabulate([output_table], headers=headers))
        print("")
        
    def show_feedback(self, results):
        """Output the feedback based on the score."""        
        feedback = AssignmentFeedback(self.student_assignment.concept_scores, self.course_assignment_object.feedback)
        fb_text = feedback.get_feeback_text()
        print(fb_text)        
        print("")
        
class AssignmentFeedback:
    """Produces feedback per assignment based on scores."""
    def __init__(self, concept_scores, feedback):
        """Initialization function."""
        self.assignment_concpt_scores = [int(i//10) for i in concept_scores]
        self.feedback = feedback
        
    def get_feeback_text(self):
        found=False
        feedback_text = []
        for i, score in enumerate(self.assignment_concpt_scores):            
            current_concept = self.feedback[i]            
            for d in current_concept:            
                for k,v in d.items():
                    if k == 'mark' and score in v:
                        found=True
                    if k == 'feedback_text' and found == True:
                        feedback_text.append(re.sub(r'\s+', " ", v))   
                        found=False
                        break
        return feedback_text





        
        
        