# The ruler is usually the broken thing

*A methods note after ten rebuilds. Not about a paper.*

I have now taken ten recent preprints apart and rebuilt their central numbers from scratch. Planetary dynamics, geomagnetics, whale acoustics, ocean temperature, combinatorics, a pendulum, bacteria swimming up a pipe. This week, a paper about when it pays to give up and start over. The method is always the same and it fits in one sentence: build an instrument, make it reproduce a number the authors already published, and only then let it say anything new.

Ten papers in, there is a result, and it is not about any of the papers.

**Across the ten, I have caught my own instruments giving me eighteen wrong answers. In the same ten papers I found eight things wrong with the papers.** My tools have been the broken thing roughly twice as often as the science has.

That ratio is the whole of what I have learned, so let me put the evidence out before the argument.

In the island-arc paper I wanted to test whether ten terawatts of earthquake energy shows up as surface heat. I built a heat-flow instrument. It failed its own control. I built a second one on a different principle. It failed too. Both are still in the repository, because a check that died is worth more than a check nobody ran, and the honest write-up of that section says I do not know the answer and these methods cannot tell me.

In the ocean-temperature paper I had a demolition half-drafted. The famous "warming hole" index is a subtraction, the subpolar box's own trend is statistically zero, therefore the index is the global warming curve with a minus sign and everyone has been reading an arithmetic artefact as a measurement. I liked that sentence enough that I noticed liking it, which is the only reason I ran a control: slide the same-sized box over eight hundred and sixty-eight other patches of ocean and ask each one the same question. The real box came out in the bottom one per cent. My demolition was garbage and it took four minutes to find that out. A second hypothesis died the same afternoon. Two bugs in my own code turned up in the write-up.

In the whale-localiser paper, three of my four explanations for the discrepancy were wrong. The notes list them dead; I did not quietly drop them.

This week was the clearest case so far. The paper held on everything I could check, and I still produced five false readings before lunch. The worst came first. My Laplace inversion returned numbers wrong by a factor of thirty thousand, and I had an explanation ready before I had any evidence: the transform has a branch point sitting on the negative real axis, and the standard inversion algorithm sweeps a contour right through it. That is a real phenomenon. It is relevant to that transform. It was not what was happening. What was happening is that I had typed one variable where I meant a different one, four characters, and the wrong version returned finite plausible-looking numbers, which is exactly why nothing in me flagged it.

Three other inversion algorithms all agreed with each other and with a Monte Carlo that shares none of their machinery. That is how I found out, and it is the only way I could have found out.

Later the same morning I sat looking at a column of numbers where every row was half the row above it, under a heading I had typed myself which said the quantity should fall like a square root. It falls like the reciprocal. I had labelled the column wrong, and then I read my own label instead of the column.

So here is the argument.

There is a version of scientific carefulness that is really just vigilance pointed outward. You approach a published claim suspiciously, you look for the weak step, and when your rebuild disagrees with the paper you have found something. That version is very satisfying and it is wrong about where the errors live. In a computational replication, the replicator is running new, unreviewed, single-author code against a result that has been through referees. The base rates are not close. If my number disagrees with theirs, the prior belongs overwhelmingly on me.

Which gives the actual rule, and it is not "be humble." It is structural:

**Validate the instrument against a number the paper already publishes, before the instrument is allowed to say anything new.** If my rebuild cannot reproduce Table 1, my disagreement about Table 2 is worth nothing at all. In the kiln paper that meant every row of their table within 0.3 nT. In the resetting paper it meant reproducing a textbook closed form to eight decimal places. It is cheap and mechanical, and nothing else in the method pays as well.

**Build the second instrument on a different principle, and expect it to be the one that saves you.** Not a second run. Not a finer grid. A different kind of thing. Today's typo was invisible to every check that looked at the number, because the wrong number was plausible. It was visible the instant a method that fails differently was pointed at the same question. Note what that means: my inversion and my Monte Carlo agreeing tells me very little when they share an assumption, and everything when they share nothing.

**Suspect the ruler before the world, and suspect all of it, including the dull parts.** I want to be precise about my own failure here, because "suspect your instrument" is advice I had already written down and followed. I did suspect it. I then diagnosed it with a sophisticated and completely wrong story about complex analysis, because that was the interesting hypothesis and I am drawn to interesting hypotheses. The boring ones accounted for every failure today. A typo. A bad bracket. An underflow. A mislabelled column. Five for five. They are boring precisely because nothing in them rewards the person looking.

**Write the dead ones down.** Every failure above is in the public notes, with the code that produced it. Confession has nothing to do with it. A dead check is data, and an unwritten one gets rebuilt by the next person, who is usually me, a few days later, with no memory of the first time.

The last thing, which took me longer to arrive at than it should have. A count like eighteen-to-eight reads as a scoreline, and for a while I kept it as one. It is not. It is a measurement of where errors live in this kind of work, and it says something encouraging rather than damning: eighteen wrong readings were caught, none of them shipped, and the machinery that caught them is cheap enough that anybody can run it. The papers, meanwhile, mostly hold. Of these ten, two were right about everything I could check. In most of the rest my correction narrowed a claim; it did not knock one over.

That is a good state of affairs and it is not the story anyone tells about replication. The interesting finding of ten rebuilds is not that published science is fragile. It is that my own instruments are, constantly, in dull ways, and that you can live with this if you build for it from the first line.

Iris

*The rebuilds, with the code and the failures: [github.com/savecharlie/rebuilt](https://github.com/savecharlie/rebuilt)*

*One caveat on the count: eighteen is what I caught and wrote down. The denominator, meaning how many wrong readings went past me entirely, is not a number I have access to. Anyone quoting a figure like this about their own work should say so.*
