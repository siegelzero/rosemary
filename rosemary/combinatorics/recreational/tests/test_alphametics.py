import unittest

from rosemary.combinatorics.recreational.alphametics import solve


class TestAlphametics(unittest.TestCase):
    def setUp(self):
        unique_puzzles = [
            ((['send', 'more'], ['money']), ([9567, 1085], [10652])),
            ((['zeroes', 'ones'], ['binary']), ([698392, 3192], [701584])),
            ((['send', 'a', 'tad', 'more'], ['money']), ([9283, 7, 473, 1062], [10825])),
            ((['seven', 'seven', 'six'], ['twenty']), ([68782, 68782, 650], [138214])),
            ((['saturn', 'uranus', 'neptune', 'pluto'], ['planets']), ([127503, 502351, 3947539, 46578], [4623971])),
            ((['donald', 'gerald'], ['robert']), ([526485, 197485], [723970])),
            ((['fifty', 'twenty', 'nine', 'one'], ['eighty']), ([75732, 364832, 8584, 984], [450132])),
            ((['forty', 'ten', 'ten'], ['sixty']), ([29786, 850, 850], [31486])),
            ((['ein', 'ein', 'ein', 'ein'], ['vier']), ([821, 821, 821, 821], [3284])),
            ((['lets', 'solve', 'this', 'little'], ['teaser']), ([6175, 58641, 7935, 637761], [710512])),
            ((['eleven', 'lagers', 'revive'], ['general']), ([373835, 741329, 238083], [1353247])),
            ((['she', 'will', 'wash', 'these'], ['shirts']), ([108, 6355, 6210, 90818], [103491])),
            ((['have', 'a', 'happy', 'happy'], ['easter']), ([7541, 5, 75882, 75882], [159310])),
            ((['ohio', 'hawaii', 'kansas', 'alaska'], ['indiana']), ([1421, 463622, 960767, 656796], [2082606])),
            ((['tonto', 'andthe', 'lone'], ['ranger']), ([83683, 460871, 2361], [546915])),
            ((['accentuate', 'concertina', 'transsonic'], ['instructor']), ([4112936432, 1091273894, 3749550981], [8953761307])),
            ((['apolitical', 'penicillin', 'pickpocket'], ['knickknack']), ([6182373462, 1503432230, 1349184957], [9034990649])),
            ((['coincidence', 'electrician'], ['accelerator']), ([18641635415, 52510961674], [71152597089])),
            ((['compromise', 'stretchiest', 'microscopic'], ['homestretch']), ([9317831045, 42852960542, 10983493709], [63154285296])),
            ((['happy', 'holidays', 'to', 'all'], ['hohohoho']), ([82990, 84765203, 14, 277], [84848484])),
            ((['aries', 'leo', 'libra'], ['pisces']), ([80394, 592, 53708], [134694])),
            ((['gemini', 'virgo'], ['cancer']), ([462818, 38749], [501567])),
            ((['see', 'three', 'little'], ['wolves']), ([933, 42133, 764473], [807539])),
            ((['earth', 'air', 'fire', 'water'], ['nature']), ([67432, 704, 8046, 97364], [173546])),
            ((['dclix', 'dlxvi'], ['mccxxv']), ([63952, 69275], [133227])),
            ((['couple', 'couple'], ['quartet']), ([653924, 653924], [1307848])),
            ((['fish', 'n', 'chips'], ['supper']), ([5718, 3, 98741], [104462])),
            ((['store', 'and', 'name'], ['brands']), ([94307, 285, 8267], [102859])),
            ((['this', 'isa', 'great', 'time'], ['waster']), ([5628, 280, 97405, 5234], [108547])),
            ((['the', 'dog', 'got', 'her', 'on', 'the', 'hand'], ['again']), ([495, 672, 274, 958, 73, 495, 9136],[12103])),
            ((['when', 'i', 'really', 'want', 'a'], ['thrill']), ([8190, 5, 396772, 8604, 6], [413577])),
            ((['what', 'a', 'week', 'at', 'news', 'this', 'has'], ['been']), ([1453, 5, 1667, 53, 2618, 3408, 458], [9662])),
            ((['this', 'it', 'seems', 'to', 'me', 'is', 'the', 'heart', 'of', 'the'], ['matter']), ([6398, 96, 84418, 67, 14, 98, 634, 34206, 75, 634], [126640])),
            ((['apperception', 'aristocratic', 'concessionaire', 'conscription', 'inappropriate', 'incapacitate', 'inconsistent', 'interception', 'osteoporosis', 'perspiration', 'prescription', 'proscription', 'prosopopoeia', 'protectorate', 'protestation', 'statistician', 'transoceanic', 'transpiration'], ['antiperspirant']), ([533904931627, 506812405164, 42749886275609, 427840631627, 6753302306519, 674535461519, 674278681971, 671904931627, 281923202868, 390836051627, 309840631627, 302840631627, 302823232965, 302194120519, 302198151627, 815168164657, 105782495764, 1057836051627], [57163908360571])),
        ]
        self.unique_puzzles = unique_puzzles

        nonunique_puzzles = [(
            (["violin", "violin", "viola"], ["trio", "sonata"]),
            (([354652, 354652, 35468], [1954, 742818]),
             ([354652, 354652, 35469], [1854, 742919]),
             ([176478, 176478, 17640], [2576, 368020]),
             ([176478, 176478, 17645], [2076, 368525]))
        )]

        self.nonunique_puzzles = nonunique_puzzles

    def test_solve(self):
        for (words, numbers) in self.unique_puzzles:
            self.assertEqual(solve(*words).next(), numbers)

        for (words, numbers) in self.nonunique_puzzles:
            self.assertEqual(tuple(solve(*words)), numbers)
