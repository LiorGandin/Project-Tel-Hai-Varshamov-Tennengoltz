# Project-Tel-Hai-Varshamov-Tennengoltz


Chandak, K. T. (2019, 6 19). Tutorial on algebraic deletion correction codes. Retrieved from ArXiv:
https://arxiv.org/pdf/1906.07887.pdf

Mahed Abroshan, R. V. (2018, 4 27). Efficient Systematic Encoding of Non-binary VT codes. Retrieved
from Arxiv: https://arxiv.org/pdf/1708.04071.pdf

Roth, R. M. (2006). Introduction to Coding Theory. Cambridge: Cambridge University Press.

Sloane, N. J. (2002, 7 22). On Single-Deletion-Correcting Codes. Retrieved from ArXiv:
https://arxiv.org/pdf/math/0207197.pdf

TENENGOLTS, G. (1984, 9). Nonbinary Codes, Correcting Single Deletion or Insertion. TRANSACTIONS
ON INFORMATION THEORY, pp. 766-770.

Univercity, H. (n.d.). Retrieved from Wyss Institute: https://wyss.harvard.edu/technology/dna-datastorage/

Varshamov, R. R., & Tenenholtz, G. M. (1965). A code for correcting a single asymmetric
error. Automatica i Telemekhanika, 26(2), 288-292.

Levenshtein, V. I. (1966, February). Binary codes capable of correcting deletions, insertions,
and reversals. In Soviet physics doklady (Vol. 10, No. 8, pp. 707-710).

all the code is written in Python3


Instructions

How to download:

git clone https://github.com/blabla


from the file Varshamov_Tennengoltz_Code tou need the class Vershamov_Tenengoltz to encode and decode:

from Varshamov_Tennengoltz_Code import Vershamov_Tenengoltz


Initializeing:

vt = Vershamov_Tenengoltz(k, q, a).


attribute of the class:

k: length of original message.
q: alphabet base for desired codeword.
a: desired syndrome value, 0 in defalt, we also aim for sum value 0 for best options. 


implemention:

codeword = vt.encode(message_to_encode)
//the codeword can be inserted 1 instance or deleted in 1 instance (1 value in the array)
decoded_message = vt.decode(codeword)
//decoded_message will be the same as message_to_encode



