# SECDED Hamming Code Implementation

This code implements a Single Error Correction, Double Error Detection (SECDED) Hamming code for any length of data.

Features

．Generates Hamming codes for error detection and correction.

．Supports data input of any length.

．Optional double error detection.

．Error injection function provided for testing the effectiveness of the Hamming code.

Usage

．Initialize the Hamming code structure using the init_hamming function.

．Encode the data using the hamming_codeword function.

．Errors can be injected into the encoded data using the inject_error function for testing purposes.

．Detect and correct errors using the hamming_get_error_info function.
