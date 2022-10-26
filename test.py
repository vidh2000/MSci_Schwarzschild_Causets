#%%
from math import nan
import numpy as np
import random as rand


text=input()

def letter_counter(word):
    
    # Store all letters in a list
    letters = []
    print("List of letters = ", letters)
    for letter in word:
        letters.append(letter)
    # Find all unique letters in the word
    characters_set = sorted(set(letters))
    # Count number of times the unique character appears in the word
    dictionary = dict()
    for character in characters_set:
        count=0
        for letter in letters:
            if character == letter:
                count +=1
        dictionary[character] = count
    return dictionary

result = letter_counter(text)
print(result)
