def solve(lhs, rhs):
    def word_to_num(word, mapping):
        return int(word.format(**mapping))

    lhs_formatted = ['{' + '}{'.join(list(word)) + '}' for word in lhs]
    rhs_formatted = ['{' + '}{'.join(list(word)) + '}' for word in rhs]

    # Build the set of all characters and the set of the first character in
    # each word.
    characters = set()
    first_chars = set()
    for word in lhs + rhs:
        characters.update(set(word))
        first_chars.add(word[0])

    num_characters = len(characters)

    # Instead of considering the given additive equation, we move all terms
    # to one side to set the equation equal to 0.

    # We define the signature of a character to be the value of the obtained by
    # letting that character have the value 1 and all other characters have
    # the value 0.
    signature = {}
    sig_map = dict(zip(characters, [0]*num_characters))

    for char in characters:
        sig_map[char] = 1
        lhs_value = sum([word_to_num(word, sig_map) for word in lhs_formatted])
        rhs_value = sum([word_to_num(word, sig_map) for word in rhs_formatted])
        signature[char] = lhs_value - rhs_value
        sig_map[char] = 0

    # As mentioned, the alphametic is solved when the total signature sum
    # is 0. We add the characters in the order that their signature
    # contributes, removing the lower place digits first.
    # i.e. we first assign values so that the 0s place of the total is 0,
    # then assign values so that the 10s place is 0, then the 100s place,
    # etc.

    # We compute the valuation of each character: This is the largest power
    # of 10 that divides the signature. We do this so it's easier to take
    # out the digits in the order described above.
    valuation = {}

    for (char, sig) in signature.iteritems():
        valuation[char] = 1
        while sig and sig % 10 == 0:
            valuation[char] *= 10
            sig //= 10

    # Finally, we sort the characters according to their valuation.
    char_list = sorted(characters, key=lambda char: valuation[char])
    mapping = {}

    def backtrack(idx, total, unused):
        if idx == num_characters:
            # If we've assigned a value to all characters, then we check
            # whether the assignment is valid, and yield it if so.
            if total == 0:
                lhs_values = [word_to_num(word, mapping) for word in lhs_formatted]
                rhs_values = [word_to_num(word, mapping) for word in rhs_formatted]
                yield lhs_values, rhs_values

        else:
            # Otherwise, we fill the character at position idx with its
            # potential values
            current_char = char_list[idx]

            for digit in unused:
                # Skip 0 if we're filling a first character
                if digit == 0 and current_char in first_chars:
                    continue

                mapping[current_char] = digit
                new_total = total + digit*signature[current_char]

                # Check if the partial solution can be extended.
                can_extend = True
                if idx < num_characters - 1:
                    next_valuation = valuation[char_list[idx + 1]]
                    if next_valuation != valuation[current_char]:
                        # If the next character only takes out
                        # higher-digits than we're currently on, this
                        # checks that all of the lower digits are zeroed
                        # out. If they're not all zeroed out, then we can't
                        # extend.
                        if new_total % next_valuation != 0:
                            can_extend = False

                if can_extend:
                    for solution in backtrack(idx + 1, new_total, unused - {digit}):
                        yield solution

            if current_char in mapping:
                del mapping[current_char]

    # We just want one solution. Remove .next() to get the generator of all
    # solutions.
    return backtrack(0, 0, set(range(10)))
