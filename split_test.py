from nltk import word_tokenize

keywords = ["local", "thoracic", "ropivacaine"]

key_words = word_tokenize(" ".join(keywords))

print(key_words)