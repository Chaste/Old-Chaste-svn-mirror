/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef _EXCEPTION_HPP_
#define _EXCEPTION_HPP_

#include <ostream>
#include <string>
#include <sstream>

#include <cfloat>
#include <climits> //For UINT_MAX etc., necessary in gcc-4.3
#include <cstdlib> //For system() etc., necessary in gcc-4.3
const unsigned UNSIGNED_UNSET=UINT_MAX;
const int INT_UNSET=INT_MAX;
const double DOUBLE_UNSET=DBL_MAX;

/**
 * Exception class.
 * All exceptions thrown by this code are currently instances of this class.
 *
 * \todo Might we want this class to inherit from STL exceptions?
 */
class Exception
{
private:
    std::string mMessage; /**< Exception message */

ACT I. Scene I.
Elsinore. A platform before the Castle.

Enter two Sentinels-[first,] Francisco, [who paces up and down
at his post; then] Bernardo, [who approaches him].

  Ber. Who's there.?
  Fran. Nay, answer me. Stand and unfold yourself.
  Ber. Long live the King!
  Fran. Bernardo?
  Ber. He.
  Fran. You come most carefully upon your hour.
  Ber. 'Tis now struck twelve. Get thee to bed, Francisco.
  Fran. For this relief much thanks. 'Tis bitter cold,
    And I am sick at heart.
  Ber. Have you had quiet guard?
  Fran. Not a mouse stirring.
  Ber. Well, good night.
    If you do meet Horatio and Marcellus,
    The rivals of my watch, bid them make haste.

                    Enter Horatio and Marcellus.  

  Fran. I think I hear them. Stand, ho! Who is there?
  Hor. Friends to this ground.
  Mar. And liegemen to the Dane.
  Fran. Give you good night.
  Mar. O, farewell, honest soldier.
    Who hath reliev'd you?
  Fran. Bernardo hath my place.
    Give you good night.                                   Exit.
  Mar. Holla, Bernardo!
  Ber. Say-
    What, is Horatio there ?
  Hor. A piece of him.
  Ber. Welcome, Horatio. Welcome, good Marcellus.
  Mar. What, has this thing appear'd again to-night?
  Ber. I have seen nothing.
  Mar. Horatio says 'tis but our fantasy,
    And will not let belief take hold of him
    Touching this dreaded sight, twice seen of us.
    Therefore I have entreated him along,  
    With us to watch the minutes of this night,
    That, if again this apparition come,
    He may approve our eyes and speak to it.
  Hor. Tush, tush, 'twill not appear.
  Ber. Sit down awhile,
    And let us once again assail your ears,
    That are so fortified against our story,
    What we two nights have seen.
  Hor. Well, sit we down,
    And let us hear Bernardo speak of this.
  Ber. Last night of all,
    When yond same star that's westward from the pole
    Had made his course t' illume that part of heaven
    Where now it burns, Marcellus and myself,
    The bell then beating one-

                        Enter Ghost.

  Mar. Peace! break thee off! Look where it comes again!
  Ber. In the same figure, like the King that's dead.  
  Mar. Thou art a scholar; speak to it, Horatio.
  Ber. Looks it not like the King? Mark it, Horatio.
  Hor. Most like. It harrows me with fear and wonder.
  Ber. It would be spoke to.
  Mar. Question it, Horatio.
  Hor. What art thou that usurp'st this time of night
    Together with that fair and warlike form
    In which the majesty of buried Denmark
    Did sometimes march? By heaven I charge thee speak!
  Mar. It is offended.
  Ber. See, it stalks away!
  Hor. Stay! Speak, speak! I charge thee speak!
                                                     Exit Ghost.
  Mar. 'Tis gone and will not answer.
  Ber. How now, Horatio? You tremble and look pale.
    Is not this something more than fantasy?
    What think you on't?
  Hor. Before my God, I might not this believe
    Without the sensible and true avouch
    Of mine own eyes.  
  Mar. Is it not like the King?
  Hor. As thou art to thyself.
    Such was the very armour he had on
    When he th' ambitious Norway combated.
    So frown'd he once when, in an angry parle,
    He smote the sledded Polacks on the ice.
    'Tis strange.
  Mar. Thus twice before, and jump at this dead hour,
    With martial stalk hath he gone by our watch.
  Hor. In what particular thought to work I know not;
    But, in the gross and scope of my opinion,
    This bodes some strange eruption to our state.
  Mar. Good now, sit down, and tell me he that knows,
    Why this same strict and most observant watch
    So nightly toils the subject of the land,
    And why such daily cast of brazen cannon
    And foreign mart for implements of war;
    Why such impress of shipwrights, whose sore task
    Does not divide the Sunday from the week.
    What might be toward, that this sweaty haste  
    Doth make the night joint-labourer with the day?
    Who is't that can inform me?
  Hor. That can I.
    At least, the whisper goes so. Our last king,
    Whose image even but now appear'd to us,
    Was, as you know, by Fortinbras of Norway,
    Thereto prick'd on by a most emulate pride,
    Dar'd to the combat; in which our valiant Hamlet
    (For so this side of our known world esteem'd him)
    Did slay this Fortinbras; who, by a seal'd compact,
    Well ratified by law and heraldry,
    Did forfeit, with his life, all those his lands
    Which he stood seiz'd of, to the conqueror;
    Against the which a moiety competent
    Was gaged by our king; which had return'd
    To the inheritance of Fortinbras,
    Had he been vanquisher, as, by the same comart
    And carriage of the article design'd,
    His fell to Hamlet. Now, sir, young Fortinbras,
    Of unimproved mettle hot and full,  
    Hath in the skirts of Norway, here and there,
    Shark'd up a list of lawless resolutes,
    For food and diet, to some enterprise
    That hath a stomach in't; which is no other,
    As it doth well appear unto our state,
    But to recover of us, by strong hand
    And terms compulsatory, those foresaid lands
    So by his father lost; and this, I take it,
    Is the main motive of our preparations,
    The source of this our watch, and the chief head
    Of this post-haste and romage in the land.
  Ber. I think it be no other but e'en so.
    Well may it sort that this portentous figure
    Comes armed through our watch, so like the King
    That was and is the question of these wars.
  Hor. A mote it is to trouble the mind's eye.
    In the most high and palmy state of Rome,
    A little ere the mightiest Julius fell,
    The graves stood tenantless, and the sheeted dead
    Did squeak and gibber in the Roman streets;  
    As stars with trains of fire, and dews of blood,
    Disasters in the sun; and the moist star
    Upon whose influence Neptune's empire stands
    Was sick almost to doomsday with eclipse.
    And even the like precurse of fierce events,
    As harbingers preceding still the fates
    And prologue to the omen coming on,
    Have heaven and earth together demonstrated
    Unto our climature and countrymen.

                      Enter Ghost again.

    But soft! behold! Lo, where it comes again!
    I'll cross it, though it blast me.- Stay illusion!
                                               Spreads his arms.
    If thou hast any sound, or use of voice,
    Speak to me.
    If there be any good thing to be done,
    That may to thee do ease, and, race to me,
    Speak to me.  
    If thou art privy to thy country's fate,
    Which happily foreknowing may avoid,
    O, speak!
    Or if thou hast uphoarded in thy life
    Extorted treasure in the womb of earth
    (For which, they say, you spirits oft walk in death),
                                                 The cock crows.
    Speak of it! Stay, and speak!- Stop it, Marcellus!
  Mar. Shall I strike at it with my partisan?
  Hor. Do, if it will not stand.
  Ber. 'Tis here!
  Hor. 'Tis here!
  Mar. 'Tis gone!
                                                     Exit Ghost.
    We do it wrong, being so majestical,
    To offer it the show of violence;
    For it is as the air, invulnerable,
    And our vain blows malicious mockery.
  Ber. It was about to speak, when the cock crew.
  Hor. And then it started, like a guilty thing  
    Upon a fearful summons. I have heard
    The cock, that is the trumpet to the morn,
    Doth with his lofty and shrill-sounding throat
    Awake the god of day; and at his warning,
    Whether in sea or fire, in earth or air,
    Th' extravagant and erring spirit hies
    To his confine; and of the truth herein
    This present object made probation.
  Mar. It faded on the crowing of the cock.
    Some say that ever, 'gainst that season comes
    Wherein our Saviour's birth is celebrated,
    The bird of dawning singeth all night long;
    And then, they say, no spirit dare stir abroad,
    The nights are wholesome, then no planets strike,
    No fairy takes, nor witch hath power to charm,
    So hallow'd and so gracious is the time.
  Hor. So have I heard and do in part believe it.
    But look, the morn, in russet mantle clad,
    Walks o'er the dew of yon high eastward hill.
    Break we our watch up; and by my advice  
    Let us impart what we have seen to-night
    Unto young Hamlet; for, upon my life,
    This spirit, dumb to us, will speak to him.
    Do you consent we shall acquaint him with it,
    As needful in our loves, fitting our duty?
    Let's do't, I pray; and I this morning know
    Where we shall find him most conveniently.           Exeunt.




Scene II.
Elsinore. A room of state in the Castle.

Flourish. [Enter Claudius, King of Denmark, Gertrude the Queen, Hamlet,
Polonius, Laertes and his sister Ophelia, [Voltemand, Cornelius,]
Lords Attendant.

  King. Though yet of Hamlet our dear brother's death
    The memory be green, and that it us befitted
    To bear our hearts in grief, and our whole kingdom
    To be contracted in one brow of woe,
    Yet so far hath discretion fought with nature
    That we with wisest sorrow think on him
    Together with remembrance of ourselves.
    Therefore our sometime sister, now our queen,
    Th' imperial jointress to this warlike state,
    Have we, as 'twere with a defeated joy,
    With an auspicious, and a dropping eye,
    With mirth in funeral, and with dirge in marriage,
    In equal scale weighing delight and dole,
    Taken to wife; nor have we herein barr'd
    Your better wisdoms, which have freely gone  
    With this affair along. For all, our thanks.
    Now follows, that you know, young Fortinbras,
    Holding a weak supposal of our worth,
    Or thinking by our late dear brother's death
    Our state to be disjoint and out of frame,
    Colleagued with this dream of his advantage,
    He hath not fail'd to pester us with message
    Importing the surrender of those lands
    Lost by his father, with all bands of law,
    To our most valiant brother. So much for him.
    Now for ourself and for this time of meeting.
    Thus much the business is: we have here writ
    To Norway, uncle of young Fortinbras,
    Who, impotent and bedrid, scarcely hears
    Of this his nephew's purpose, to suppress
    His further gait herein, in that the levies,
    The lists, and full proportions are all made
    Out of his subject; and we here dispatch
    You, good Cornelius, and you, Voltemand,
    For bearers of this greeting to old Norway,  
    Giving to you no further personal power
    To business with the King, more than the scope
    Of these dilated articles allow.            [Gives a paper.]
    Farewell, and let your haste commend your duty.
  Cor., Volt. In that, and all things, will we show our duty.
  King. We doubt it nothing. Heartily farewell.
                                 Exeunt Voltemand and Cornelius.
    And now, Laertes, what's the news with you?
    You told us of some suit. What is't, Laertes?
    You cannot speak of reason to the Dane
    And lose your voice. What wouldst thou beg, Laertes,
    That shall not be my offer, not thy asking?
    The head is not more native to the heart,
    The hand more instrumental to the mouth,
    Than is the throne of Denmark to thy father.
    What wouldst thou have, Laertes?
  Laer. My dread lord,
    Your leave and favour to return to France;
    From whence though willingly I came to Denmark
    To show my duty in your coronation,  
    Yet now I must confess, that duty done,
    My thoughts and wishes bend again toward France
    And bow them to your gracious leave and pardon.
  King. Have you your father's leave? What says Polonius?
  Pol. He hath, my lord, wrung from me my slow leave
    By laboursome petition, and at last
    Upon his will I seal'd my hard consent.
    I do beseech you give him leave to go.
  King. Take thy fair hour, Laertes. Time be thine,
    And thy best graces spend it at thy will!
    But now, my cousin Hamlet, and my son-
  Ham. [aside] A little more than kin, and less than kind!
  King. How is it that the clouds still hang on you?
  Ham. Not so, my lord. I am too much i' th' sun.
  Queen. Good Hamlet, cast thy nighted colour off,
    And let thine eye look like a friend on Denmark.
    Do not for ever with thy vailed lids
    Seek for thy noble father in the dust.
    Thou know'st 'tis common. All that lives must die,
    Passing through nature to eternity.  
  Ham. Ay, madam, it is common.
  Queen. If it be,
    Why seems it so particular with thee?
  Ham. Seems, madam, Nay, it is. I know not 'seems.'
    'Tis not alone my inky cloak, good mother,
    Nor customary suits of solemn black,
    Nor windy suspiration of forc'd breath,
    No, nor the fruitful river in the eye,
    Nor the dejected havior of the visage,
    Together with all forms, moods, shapes of grief,
    'That can denote me truly. These indeed seem,
    For they are actions that a man might play;
    But I have that within which passeth show-
    These but the trappings and the suits of woe.
  King. 'Tis sweet and commendable in your nature, Hamlet,
    To give these mourning duties to your father;
    But you must know, your father lost a father;
    That father lost, lost his, and the survivor bound
    In filial obligation for some term
    To do obsequious sorrow. But to persever  
    In obstinate condolement is a course
    Of impious stubbornness. 'Tis unmanly grief;
    It shows a will most incorrect to heaven,
    A heart unfortified, a mind impatient,
    An understanding simple and unschool'd;
    For what we know must be, and is as common
    As any the most vulgar thing to sense,
    Why should we in our peevish opposition
    Take it to heart? Fie! 'tis a fault to heaven,
    A fault against the dead, a fault to nature,
    To reason most absurd, whose common theme
    Is death of fathers, and who still hath cried,
    From the first corse till he that died to-day,
    'This must be so.' We pray you throw to earth
    This unprevailing woe, and think of us
    As of a father; for let the world take note
    You are the most immediate to our throne,
    And with no less nobility of love
    Than that which dearest father bears his son
    Do I impart toward you. For your intent  
    In going back to school in Wittenberg,
    It is most retrograde to our desire;
    And we beseech you, bend you to remain
    Here in the cheer and comfort of our eye,
    Our chiefest courtier, cousin, and our son.
  Queen. Let not thy mother lose her prayers, Hamlet.
    I pray thee stay with us, go not to Wittenberg.
  Ham. I shall in all my best obey you, madam.
  King. Why, 'tis a loving and a fair reply.
    Be as ourself in Denmark. Madam, come.
    This gentle and unforc'd accord of Hamlet
    Sits smiling to my heart; in grace whereof,
    No jocund health that Denmark drinks to-day
    But the great cannon to the clouds shall tell,
    And the King's rouse the heaven shall bruit again,
    Respeaking earthly thunder. Come away.
                                Flourish. Exeunt all but Hamlet.
  Ham. O that this too too solid flesh would melt,
    Thaw, and resolve itself into a dew!
    Or that the Everlasting had not fix'd  
    His canon 'gainst self-slaughter! O God! God!
    How weary, stale, flat, and unprofitable
    Seem to me all the uses of this world!
    Fie on't! ah, fie! 'Tis an unweeded garden
    That grows to seed; things rank and gross in nature
    Possess it merely. That it should come to this!
    But two months dead! Nay, not so much, not two.
    So excellent a king, that was to this
    Hyperion to a satyr; so loving to my mother
    That he might not beteem the winds of heaven
    Visit her face too roughly. Heaven and earth!
    Must I remember? Why, she would hang on him
    As if increase of appetite had grown
    By what it fed on; and yet, within a month-
    Let me not think on't! Frailty, thy name is woman!-
    A little month, or ere those shoes were old
    With which she followed my poor father's body
    Like Niobe, all tears- why she, even she
    (O God! a beast that wants discourse of reason
    Would have mourn'd longer) married with my uncle;  
    My father's brother, but no more like my father
    Than I to Hercules. Within a month,
    Ere yet the salt of most unrighteous tears
    Had left the flushing in her galled eyes,
    She married. O, most wicked speed, to post
    With such dexterity to incestuous sheets!
    It is not, nor it cannot come to good.
    But break my heart, for I must hold my tongue!

          Enter Horatio, Marcellus, and Bernardo.

  Hor. Hail to your lordship!
  Ham. I am glad to see you well.
    Horatio!- or I do forget myself.
  Hor. The same, my lord, and your poor servant ever.
  Ham. Sir, my good friend- I'll change that name with you.
    And what make you from Wittenberg, Horatio?
    Marcellus?
  Mar. My good lord!
  Ham. I am very glad to see you.- [To Bernardo] Good even, sir.-  
    But what, in faith, make you from Wittenberg?
  Hor. A truant disposition, good my lord.
  Ham. I would not hear your enemy say so,
    Nor shall you do my ear that violence
    To make it truster of your own report
    Against yourself. I know you are no truant.
    But what is your affair in Elsinore?
    We'll teach you to drink deep ere you depart.
  Hor. My lord, I came to see your father's funeral.
  Ham. I prithee do not mock me, fellow student.
    I think it was to see my mother's wedding.
  Hor. Indeed, my lord, it followed hard upon.
  Ham. Thrift, thrift, Horatio! The funeral bak'd meats
    Did coldly furnish forth the marriage tables.
    Would I had met my dearest foe in heaven
    Or ever I had seen that day, Horatio!
    My father- methinks I see my father.
  Hor. O, where, my lord?
  Ham. In my mind's eye, Horatio.
  Hor. I saw him once. He was a goodly king.  
  Ham. He was a man, take him for all in all.
    I shall not look upon his like again.
  Hor. My lord, I think I saw him yesternight.
  Ham. Saw? who?
  Hor. My lord, the King your father.
  Ham. The King my father?
  Hor. Season your admiration for a while
    With an attent ear, till I may deliver
    Upon the witness of these gentlemen,
    This marvel to you.
  Ham. For God's love let me hear!
  Hor. Two nights together had these gentlemen
    (Marcellus and Bernardo) on their watch
    In the dead vast and middle of the night
    Been thus encount'red. A figure like your father,
    Armed at point exactly, cap-a-pe,
    Appears before them and with solemn march
    Goes slow and stately by them. Thrice he walk'd
    By their oppress'd and fear-surprised eyes,
    Within his truncheon's length; whilst they distill'd  
    Almost to jelly with the act of fear,
    Stand dumb and speak not to him. This to me
    In dreadful secrecy impart they did,
    And I with them the third night kept the watch;
    Where, as they had deliver'd, both in time,
    Form of the thing, each word made true and good,
    The apparition comes. I knew your father.
    These hands are not more like.
  Ham. But where was this?
  Mar. My lord, upon the platform where we watch'd.
  Ham. Did you not speak to it?
  Hor. My lord, I did;
    But answer made it none. Yet once methought
    It lifted up it head and did address
    Itself to motion, like as it would speak;
    But even then the morning cock crew loud,
    And at the sound it shrunk in haste away
    And vanish'd from our sight.
  Ham. 'Tis very strange.
  Hor. As I do live, my honour'd lord, 'tis true;  
    And we did think it writ down in our duty
    To let you know of it.
  Ham. Indeed, indeed, sirs. But this troubles me.
    Hold you the watch to-night?
  Both [Mar. and Ber.] We do, my lord.
  Ham. Arm'd, say you?
  Both. Arm'd, my lord.
  Ham. From top to toe?
  Both. My lord, from head to foot.
  Ham. Then saw you not his face?
  Hor. O, yes, my lord! He wore his beaver up.
  Ham. What, look'd he frowningly.
  Hor. A countenance more in sorrow than in anger.
  Ham. Pale or red?
  Hor. Nay, very pale.
  Ham. And fix'd his eyes upon you?
  Hor. Most constantly.
  Ham. I would I had been there.
  Hor. It would have much amaz'd you.
  Ham. Very like, very like. Stay'd it long?  
  Hor. While one with moderate haste might tell a hundred.
  Both. Longer, longer.
  Hor. Not when I saw't.
  Ham. His beard was grizzled- no?
  Hor. It was, as I have seen it in his life,
    A sable silver'd.
  Ham. I will watch to-night.
    Perchance 'twill walk again.
  Hor. I warr'nt it will.
  Ham. If it assume my noble father's person,
    I'll speak to it, though hell itself should gape
    And bid me hold my peace. I pray you all,
    If you have hitherto conceal'd this sight,
    Let it be tenable in your silence still;
    And whatsoever else shall hap to-night,
    Give it an understanding but no tongue.
    I will requite your loves. So, fare you well.
    Upon the platform, 'twixt eleven and twelve,
    I'll visit you.
  All. Our duty to your honour.  
  Ham. Your loves, as mine to you. Farewell.
                                        Exeunt [all but Hamlet].
    My father's spirit- in arms? All is not well.
    I doubt some foul play. Would the night were come!
    Till then sit still, my soul. Foul deeds will rise,
    Though all the earth o'erwhelm them, to men's eyes.
Exit.




Scene III.
Elsinore. A room in the house of Polonius.

Enter Laertes and Ophelia.

  Laer. My necessaries are embark'd. Farewell.
    And, sister, as the winds give benefit
    And convoy is assistant, do not sleep,
    But let me hear from you.
  Oph. Do you doubt that?
  Laer. For Hamlet, and the trifling of his favour,
    Hold it a fashion, and a toy in blood;
    A violet in the youth of primy nature,
    Forward, not permanent- sweet, not lasting;
    The perfume and suppliance of a minute;
    No more.
  Oph. No more but so?
  Laer. Think it no more.
    For nature crescent does not grow alone
    In thews and bulk; but as this temple waxes,
    The inward service of the mind and soul
    Grows wide withal. Perhaps he loves you now,  
    And now no soil nor cautel doth besmirch
    The virtue of his will; but you must fear,
    His greatness weigh'd, his will is not his own;
    For he himself is subject to his birth.
    He may not, as unvalued persons do,
    Carve for himself, for on his choice depends
    The safety and health of this whole state,
    And therefore must his choice be circumscrib'd
    Unto the voice and yielding of that body
    Whereof he is the head. Then if he says he loves you,
    It fits your wisdom so far to believe it
    As he in his particular act and place
    May give his saying deed; which is no further
    Than the main voice of Denmark goes withal.
    Then weigh what loss your honour may sustain
    If with too credent ear you list his songs,
    Or lose your heart, or your chaste treasure open
    To his unmast'red importunity.
    Fear it, Ophelia, fear it, my dear sister,
    And keep you in the rear of your affection,  
    Out of the shot and danger of desire.
    The chariest maid is prodigal enough
    If she unmask her beauty to the moon.
    Virtue itself scopes not calumnious strokes.
    The canker galls the infants of the spring
    Too oft before their buttons be disclos'd,
    And in the morn and liquid dew of youth
    Contagious blastments are most imminent.
    Be wary then; best safety lies in fear.
    Youth to itself rebels, though none else near.
  Oph. I shall th' effect of this good lesson keep
    As watchman to my heart. But, good my brother,
    Do not as some ungracious pastors do,
    Show me the steep and thorny way to heaven,
    Whiles, like a puff'd and reckless libertine,
    Himself the primrose path of dalliance treads
    And recks not his own rede.
  Laer. O, fear me not!

                       Enter Polonius.  

    I stay too long. But here my father comes.
    A double blessing is a double grace;
    Occasion smiles upon a second leave.
  Pol. Yet here, Laertes? Aboard, aboard, for shame!
    The wind sits in the shoulder of your sail,
    And you are stay'd for. There- my blessing with thee!
    And these few precepts in thy memory
    Look thou character. Give thy thoughts no tongue,
    Nor any unproportion'd thought his act.
    Be thou familiar, but by no means vulgar:
    Those friends thou hast, and their adoption tried,
    Grapple them unto thy soul with hoops of steel;
    But do not dull thy palm with entertainment
    Of each new-hatch'd, unfledg'd comrade. Beware
    Of entrance to a quarrel; but being in,
    Bear't that th' opposed may beware of thee.
    Give every man thine ear, but few thy voice;
    Take each man's censure, but reserve thy judgment.
    Costly thy habit as thy purse can buy,  
    But not express'd in fancy; rich, not gaudy;
    For the apparel oft proclaims the man,
    And they in France of the best rank and station
    Are most select and generous, chief in that.
    Neither a borrower nor a lender be;
    For loan oft loses both itself and friend,
    And borrowing dulls the edge of husbandry.
    This above all- to thine own self be true,
    And it must follow, as the night the day,
    Thou canst not then be false to any man.
    Farewell. My blessing season this in thee!
  Laer. Most humbly do I take my leave, my lord.
  Pol. The time invites you. Go, your servants tend.
  Laer. Farewell, Ophelia, and remember well
    What I have said to you.
  Oph. 'Tis in my memory lock'd,
    And you yourself shall keep the key of it.
  Laer. Farewell.                                          Exit.
  Pol. What is't, Ophelia, he hath said to you?
  Oph. So please you, something touching the Lord Hamlet.  
  Pol. Marry, well bethought!
    'Tis told me he hath very oft of late
    Given private time to you, and you yourself
    Have of your audience been most free and bounteous.
    If it be so- as so 'tis put on me,
    And that in way of caution- I must tell you
    You do not understand yourself so clearly
    As it behooves my daughter and your honour.
    What is between you? Give me up the truth.
  Oph. He hath, my lord, of late made many tenders
    Of his affection to me.
  Pol. Affection? Pooh! You speak like a green girl,
    Unsifted in such perilous circumstance.
    Do you believe his tenders, as you call them?
  Oph. I do not know, my lord, what I should think,
  Pol. Marry, I will teach you! Think yourself a baby
    That you have ta'en these tenders for true pay,
    Which are not sterling. Tender yourself more dearly,
    Or (not to crack the wind of the poor phrase,
    Running it thus) you'll tender me a fool.  
  Oph. My lord, he hath importun'd me with love
    In honourable fashion.
  Pol. Ay, fashion you may call it. Go to, go to!
  Oph. And hath given countenance to his speech, my lord,
    With almost all the holy vows of heaven.
  Pol. Ay, springes to catch woodcocks! I do know,
    When the blood burns, how prodigal the soul
    Lends the tongue vows. These blazes, daughter,
    Giving more light than heat, extinct in both
    Even in their promise, as it is a-making,
    You must not take for fire. From this time
    Be something scanter of your maiden presence.
    Set your entreatments at a higher rate
    Than a command to parley. For Lord Hamlet,
    Believe so much in him, that he is young,
    And with a larger tether may he walk
    Than may be given you. In few, Ophelia,
    Do not believe his vows; for they are brokers,
    Not of that dye which their investments show,
    But mere implorators of unholy suits,  
    Breathing like sanctified and pious bawds,
    The better to beguile. This is for all:
    I would not, in plain terms, from this time forth
    Have you so slander any moment leisure
    As to give words or talk with the Lord Hamlet.
    Look to't, I charge you. Come your ways.
  Oph. I shall obey, my lord.
                                                         Exeunt.




Scene IV.
Elsinore. The platform before the Castle.

Enter Hamlet, Horatio, and Marcellus.

  Ham. The air bites shrewdly; it is very cold.
  Hor. It is a nipping and an eager air.
  Ham. What hour now?
  Hor. I think it lacks of twelve.
  Mar. No, it is struck.
  Hor. Indeed? I heard it not. It then draws near the season
    Wherein the spirit held his wont to walk.
                   A flourish of trumpets, and two pieces go off.
    What does this mean, my lord?
  Ham. The King doth wake to-night and takes his rouse,
    Keeps wassail, and the swagg'ring upspring reels,
    And, as he drains his draughts of Rhenish down,
    The kettledrum and trumpet thus bray out
    The triumph of his pledge.
  Hor. Is it a custom?
  Ham. Ay, marry, is't;
    But to my mind, though I am native here  
    And to the manner born, it is a custom
    More honour'd in the breach than the observance.
    This heavy-headed revel east and west
    Makes us traduc'd and tax'd of other nations;
    They clip us drunkards and with swinish phrase
    Soil our addition; and indeed it takes
    From our achievements, though perform'd at height,
    The pith and marrow of our attribute.
    So oft it chances in particular men
    That, for some vicious mole of nature in them,
    As in their birth,- wherein they are not guilty,
    Since nature cannot choose his origin,-
    By the o'ergrowth of some complexion,
    Oft breaking down the pales and forts of reason,
    Or by some habit that too much o'erleavens
    The form of plausive manners, that these men
    Carrying, I say, the stamp of one defect,
    Being nature's livery, or fortune's star,
    Their virtues else- be they as pure as grace,
    As infinite as man may undergo-  
    Shall in the general censure take corruption
    From that particular fault. The dram of e'il
    Doth all the noble substance often dout To his own scandal.

                         Enter Ghost.

  Hor. Look, my lord, it comes!
  Ham. Angels and ministers of grace defend us!
    Be thou a spirit of health or goblin damn'd,
    Bring with thee airs from heaven or blasts from hell,
    Be thy intents wicked or charitable,
    Thou com'st in such a questionable shape
    That I will speak to thee. I'll call thee Hamlet,
    King, father, royal Dane. O, answer me?
    Let me not burst in ignorance, but tell
    Why thy canoniz'd bones, hearsed in death,
    Have burst their cerements; why the sepulchre
    Wherein we saw thee quietly inurn'd,
    Hath op'd his ponderous and marble jaws
    To cast thee up again. What may this mean  
    That thou, dead corse, again in complete steel,
    Revisits thus the glimpses of the moon,
    Making night hideous, and we fools of nature
    So horridly to shake our disposition
    With thoughts beyond the reaches of our souls?
    Say, why is this? wherefore? What should we do?
                                           Ghost beckons Hamlet.
  Hor. It beckons you to go away with it,
    As if it some impartment did desire
    To you alone.
  Mar. Look with what courteous action
    It waves you to a more removed ground.
    But do not go with it!
  Hor. No, by no means!
  Ham. It will not speak. Then will I follow it.
  Hor. Do not, my lord!
  Ham. Why, what should be the fear?
    I do not set my life at a pin's fee;
    And for my soul, what can it do to that,
    Being a thing immortal as itself?  
    It waves me forth again. I'll follow it.
  Hor. What if it tempt you toward the flood, my lord,
    Or to the dreadful summit of the cliff
    That beetles o'er his base into the sea,
    And there assume some other, horrible form
    Which might deprive your sovereignty of reason
    And draw you into madness? Think of it.
    The very place puts toys of desperation,
    Without more motive, into every brain
    That looks so many fadoms to the sea
    And hears it roar beneath.
  Ham. It waves me still.
    Go on. I'll follow thee.
  Mar. You shall not go, my lord.
  Ham. Hold off your hands!
  Hor. Be rul'd. You shall not go.
  Ham. My fate cries out
    And makes each petty artire in this body
    As hardy as the Nemean lion's nerve.
                                                [Ghost beckons.]  
    Still am I call'd. Unhand me, gentlemen.
    By heaven, I'll make a ghost of him that lets me!-
    I say, away!- Go on. I'll follow thee.
                                        Exeunt Ghost and Hamlet.
  Hor. He waxes desperate with imagination.
  Mar. Let's follow. 'Tis not fit thus to obey him.
  Hor. Have after. To what issue wail this come?
  Mar. Something is rotten in the state of Denmark.
  Hor. Heaven will direct it.
  Mar. Nay, let's follow him.
                                                         Exeunt.




Scene V.
Elsinore. The Castle. Another part of the fortifications.

Enter Ghost and Hamlet.

  Ham. Whither wilt thou lead me? Speak! I'll go no further.
  Ghost. Mark me.
  Ham. I will.
  Ghost. My hour is almost come,
    When I to sulph'rous and tormenting flames
    Must render up myself.
  Ham. Alas, poor ghost!
  Ghost. Pity me not, but lend thy serious hearing
    To what I shall unfold.
  Ham. Speak. I am bound to hear.
  Ghost. So art thou to revenge, when thou shalt hear.
  Ham. What?
  Ghost. I am thy father's spirit,
    Doom'd for a certain term to walk the night,
    And for the day confin'd to fast in fires,
    Till the foul crimes done in my days of nature
    Are burnt and purg'd away. But that I am forbid  
    To tell the secrets of my prison house,
    I could a tale unfold whose lightest word
    Would harrow up thy soul, freeze thy young blood,
    Make thy two eyes, like stars, start from their spheres,
    Thy knotted and combined locks to part,
    And each particular hair to stand an end
    Like quills upon the fretful porpentine.
    But this eternal blazon must not be
    To ears of flesh and blood. List, list, O, list!
    If thou didst ever thy dear father love-
  Ham. O God!
  Ghost. Revenge his foul and most unnatural murther.
  Ham. Murther?
  Ghost. Murther most foul, as in the best it is;
    But this most foul, strange, and unnatural.
  Ham. Haste me to know't, that I, with wings as swift
    As meditation or the thoughts of love,
    May sweep to my revenge.
  Ghost. I find thee apt;
    And duller shouldst thou be than the fat weed  
    That rots itself in ease on Lethe wharf,
    Wouldst thou not stir in this. Now, Hamlet, hear.
    'Tis given out that, sleeping in my orchard,
    A serpent stung me. So the whole ear of Denmark
    Is by a forged process of my death
    Rankly abus'd. But know, thou noble youth,
    The serpent that did sting thy father's life
    Now wears his crown.
  Ham. O my prophetic soul!
    My uncle?
  Ghost. Ay, that incestuous, that adulterate beast,
    With witchcraft of his wit, with traitorous gifts-
    O wicked wit and gifts, that have the power
    So to seduce!- won to his shameful lust
    The will of my most seeming-virtuous queen.
    O Hamlet, what a falling-off was there,
    From me, whose love was of that dignity
    That it went hand in hand even with the vow
    I made to her in marriage, and to decline
    Upon a wretch whose natural gifts were poor  
    To those of mine!
    But virtue, as it never will be mov'd,
    Though lewdness court it in a shape of heaven,
    So lust, though to a radiant angel link'd,
    Will sate itself in a celestial bed
    And prey on garbage.
    But soft! methinks I scent the morning air.
    Brief let me be. Sleeping within my orchard,
    My custom always of the afternoon,
    Upon my secure hour thy uncle stole,
    With juice of cursed hebona in a vial,
    And in the porches of my ears did pour
    The leperous distilment; whose effect
    Holds such an enmity with blood of man
    That swift as quicksilverr it courses through
    The natural gates and alleys of the body,
    And with a sudden vigour it doth posset
    And curd, like eager droppings into milk,
    The thin and wholesome blood. So did it mine;
    And a most instant tetter bark'd about,  
    Most lazar-like, with vile and loathsome crust
    All my smooth body.
    Thus was I, sleeping, by a brother's hand
    Of life, of crown, of queen, at once dispatch'd;
    Cut off even in the blossoms of my sin,
    Unhous'led, disappointed, unanel'd,
    No reckoning made, but sent to my account
    With all my imperfections on my head.
  Ham. O, horrible! O, horrible! most horrible!
  Ghost. If thou hast nature in thee, bear it not.
    Let not the royal bed of Denmark be
    A couch for luxury and damned incest.
    But, howsoever thou pursuest this act,
    Taint not thy mind, nor let thy soul contrive
    Against thy mother aught. Leave her to heaven,
    And to those thorns that in her bosom lodge
    To prick and sting her. Fare thee well at once.
    The glowworm shows the matin to be near
    And gins to pale his uneffectual fire.
    Adieu, adieu, adieu! Remember me.                      Exit.  
  Ham. O all you host of heaven! O earth! What else?
    And shall I couple hell? Hold, hold, my heart!
    And you, my sinews, grow not instant old,
    But bear me stiffly up. Remember thee?
    Ay, thou poor ghost, while memory holds a seat
    In this distracted globe. Remember thee?
    Yea, from the table of my memory
    I'll wipe away all trivial fond records,
    All saws of books, all forms, all pressures past
    That youth and observation copied there,
    And thy commandment all alone shall live
    Within the book and volume of my brain,
    Unmix'd with baser matter. Yes, by heaven!
    O most pernicious woman!
    O villain, villain, smiling, damned villain!
    My tables! Meet it is I set it down
    That one may smile, and smile, and be a villain;
    At least I am sure it may be so in Denmark.        [Writes.]
    So, uncle, there you are. Now to my word:
    It is 'Adieu, adieu! Remember me.'  
    I have sworn't.
  Hor. (within) My lord, my lord!

                   Enter Horatio and Marcellus.

  Mar. Lord Hamlet!
  Hor. Heaven secure him!
  Ham. So be it!
  Mar. Illo, ho, ho, my lord!
  Ham. Hillo, ho, ho, boy! Come, bird, come.
  Mar. How is't, my noble lord?
  Hor. What news, my lord?
  Mar. O, wonderful!
  Hor. Good my lord, tell it.
  Ham. No, you will reveal it.
  Hor. Not I, my lord, by heaven!
  Mar. Nor I, my lord.
  Ham. How say you then? Would heart of man once think it?
    But you'll be secret?
  Both. Ay, by heaven, my lord.  
  Ham. There's neer a villain dwelling in all Denmark
    But he's an arrant knave.
  Hor. There needs no ghost, my lord, come from the grave
    To tell us this.
  Ham. Why, right! You are in the right!
    And so, without more circumstance at all,
    I hold it fit that we shake hands and part;
    You, as your business and desires shall point you,
    For every man hath business and desire,
    Such as it is; and for my own poor part,
    Look you, I'll go pray.
  Hor. These are but wild and whirling words, my lord.
  Ham. I am sorry they offend you, heartily;
    Yes, faith, heartily.
  Hor. There's no offence, my lord.
  Ham. Yes, by Saint Patrick, but there is, Horatio,
    And much offence too. Touching this vision here,
    It is an honest ghost, that let me tell you.
    For your desire to know what is between us,
    O'ermaster't as you may. And now, good friends,  
    As you are friends, scholars, and soldiers,
    Give me one poor request.
  Hor. What is't, my lord? We will.
  Ham. Never make known what you have seen to-night.
  Both. My lord, we will not.
  Ham. Nay, but swear't.
  Hor. In faith,
    My lord, not I.
  Mar. Nor I, my lord- in faith.
  Ham. Upon my sword.
  Mar. We have sworn, my lord, already.
  Ham. Indeed, upon my sword, indeed.

                 Ghost cries under the stage.

  Ghost. Swear.
  Ham. Aha boy, say'st thou so? Art thou there, truepenny?
    Come on! You hear this fellow in the cellarage.
    Consent to swear.
  Hor. Propose the oath, my lord.  
  Ham. Never to speak of this that you have seen.
    Swear by my sword.
  Ghost. [beneath] Swear.
  Ham. Hic et ubique? Then we'll shift our ground.
    Come hither, gentlemen,
    And lay your hands again upon my sword.
    Never to speak of this that you have heard:
    Swear by my sword.
  Ghost. [beneath] Swear by his sword.
  Ham. Well said, old mole! Canst work i' th' earth so fast?
    A worthy pioner! Once more remove, good friends."
  Hor. O day and night, but this is wondrous strange!
  Ham. And therefore as a stranger give it welcome.
    There are more things in heaven and earth, Horatio,
    Than are dreamt of in your philosophy.
    But come!
    Here, as before, never, so help you mercy,
    How strange or odd soe'er I bear myself
    (As I perchance hereafter shall think meet
    To put an antic disposition on),  
    That you, at such times seeing me, never shall,
    With arms encumb'red thus, or this head-shake,
    Or by pronouncing of some doubtful phrase,
    As 'Well, well, we know,' or 'We could, an if we would,'
    Or 'If we list to speak,' or 'There be, an if they might,'
    Or such ambiguous giving out, to note
    That you know aught of me- this is not to do,
    So grace and mercy at your most need help you,
    Swear.
  Ghost. [beneath] Swear.
                                                   [They swear.]
  Ham. Rest, rest, perturbed spirit! So, gentlemen,
    With all my love I do commend me to you;
    And what so poor a man as Hamlet is
    May do t' express his love and friending to you,
    God willing, shall not lack. Let us go in together;
    And still your fingers on your lips, I pray.
    The time is out of joint. O cursed spite
    That ever I was born to set it right!
    Nay, come, let's go together.  
                                                         Exeunt.




Act II. Scene I.
Elsinore. A room in the house of Polonius.

Enter Polonius and Reynaldo.

  Pol. Give him this money and these notes, Reynaldo.
  Rey. I will, my lord.
  Pol. You shall do marvell's wisely, good Reynaldo,
    Before You visit him, to make inquire
    Of his behaviour.
  Rey. My lord, I did intend it.
  Pol. Marry, well said, very well said. Look you, sir,
    Enquire me first what Danskers are in Paris;
    And how, and who, what means, and where they keep,
    What company, at what expense; and finding
    By this encompassment and drift of question
    That they do know my son, come you more nearer
    Than your particular demands will touch it.
    Take you, as 'twere, some distant knowledge of him;
    As thus, 'I know his father and his friends,
    And in part him.' Do you mark this, Reynaldo?
  Rey. Ay, very well, my lord.  
  Pol. 'And in part him, but,' you may say, 'not well.
    But if't be he I mean, he's very wild
    Addicted so and so'; and there put on him
    What forgeries you please; marry, none so rank
    As may dishonour him- take heed of that;
    But, sir, such wanton, wild, and usual slips
    As are companions noted and most known
    To youth and liberty.
  Rey. As gaming, my lord.
  Pol. Ay, or drinking, fencing, swearing, quarrelling,
    Drabbing. You may go so far.
  Rey. My lord, that would dishonour him.
  Pol. Faith, no, as you may season it in the charge.
    You must not put another scandal on him,
    That he is open to incontinency.
    That's not my meaning. But breathe his faults so quaintly
    That they may seem the taints of liberty,
    The flash and outbreak of a fiery mind,
    A savageness in unreclaimed blood,
    Of general assault.  
  Rey. But, my good lord-
  Pol. Wherefore should you do this?
  Rey. Ay, my lord,
    I would know that.
  Pol. Marry, sir, here's my drift,
    And I believe it is a fetch of warrant.
    You laying these slight sullies on my son
    As 'twere a thing a little soil'd i' th' working,
    Mark you,
    Your party in converse, him you would sound,
    Having ever seen in the prenominate crimes
    The youth you breathe of guilty, be assur'd
    He closes with you in this consequence:
    'Good sir,' or so, or 'friend,' or 'gentleman'-
    According to the phrase or the addition
    Of man and country-
  Rey. Very good, my lord.
  Pol. And then, sir, does 'a this- 'a does- What was I about to say?
    By the mass, I was about to say something! Where did I leave?
  Rey. At 'closes in the consequence,' at 'friend or so,' and  
    gentleman.'
  Pol. At 'closes in the consequence'- Ay, marry!
    He closes thus: 'I know the gentleman.
    I saw him yesterday, or t'other day,
    Or then, or then, with such or such; and, as you say,
    There was 'a gaming; there o'ertook in's rouse;
    There falling out at tennis'; or perchance,
    'I saw him enter such a house of sale,'
    Videlicet, a brothel, or so forth.
    See you now-
    Your bait of falsehood takes this carp of truth;
    And thus do we of wisdom and of reach,
    With windlasses and with assays of bias,
    By indirections find directions out.
    So, by my former lecture and advice,
    Shall you my son. You have me, have you not
  Rey. My lord, I have.
  Pol. God b' wi' ye, fare ye well!
  Rey. Good my lord!                                    [Going.]
  Pol. Observe his inclination in yourself.  
  Rey. I shall, my lord.
  Pol. And let him ply his music.
  Rey. Well, my lord.
  Pol. Farewell!
                                                  Exit Reynaldo.

                       Enter Ophelia.

    How now, Ophelia? What's the matter?
  Oph. O my lord, my lord, I have been so affrighted!
  Pol. With what, i' th' name of God I
  Oph. My lord, as I was sewing in my closet,
    Lord Hamlet, with his doublet all unbrac'd,
    No hat upon his head, his stockings foul'd,
    Ungart'red, and down-gyved to his ankle;
    Pale as his shirt, his knees knocking each other,
    And with a look so piteous in purport
    As if he had been loosed out of hell
    To speak of horrors- he comes before me.
  Pol. Mad for thy love?  
  Oph. My lord, I do not know,
    But truly I do fear it.
  Pol. What said he?
  Oph. He took me by the wrist and held me hard;
    Then goes he to the length of all his arm,
    And, with his other hand thus o'er his brow,
    He falls to such perusal of my face
    As he would draw it. Long stay'd he so.
    At last, a little shaking of mine arm,
    And thrice his head thus waving up and down,
    He rais'd a sigh so piteous and profound
    As it did seem to shatter all his bulk
    And end his being. That done, he lets me go,
    And with his head over his shoulder turn'd
    He seem'd to find his way without his eyes,
    For out o' doors he went without their help
    And to the last bended their light on me.
  Pol. Come, go with me. I will go seek the King.
    This is the very ecstasy of love,
    Whose violent property fordoes itself  
    And leads the will to desperate undertakings
    As oft as any passion under heaven
    That does afflict our natures. I am sorry.
    What, have you given him any hard words of late?
  Oph. No, my good lord; but, as you did command,
    I did repel his letters and denied
    His access to me.
  Pol. That hath made him mad.
    I am sorry that with better heed and judgment
    I had not quoted him. I fear'd he did but trifle
    And meant to wrack thee; but beshrew my jealousy!
    By heaven, it is as proper to our age
    To cast beyond ourselves in our opinions
    As it is common for the younger sort
    To lack discretion. Come, go we to the King.
    This must be known; which, being kept close, might move
    More grief to hide than hate to utter love.
    Come.
                                                         Exeunt.

Scene II.
Elsinore. A room in the Castle.

Flourish. [Enter King and Queen, Rosencrantz and Guildenstern, cum aliis.

  King. Welcome, dear Rosencrantz and Guildenstern.
    Moreover that we much did long to see you,
    The need we have to use you did provoke
    Our hasty sending. Something have you heard
    Of Hamlet's transformation. So I call it,
    Sith nor th' exterior nor the inward man
    Resembles that it was. What it should be,
    More than his father's death, that thus hath put him
    So much from th' understanding of himself,
    I cannot dream of. I entreat you both
    That, being of so young clays brought up with him,
    And since so neighbour'd to his youth and haviour,
    That you vouchsafe your rest here in our court
    Some little time; so by your companies
    To draw him on to pleasures, and to gather
    So much as from occasion you may glean,  
    Whether aught to us unknown afflicts him thus
    That, open'd, lies within our remedy.
  Queen. Good gentlemen, he hath much talk'd of you,
    And sure I am two men there are not living
    To whom he more adheres. If it will please you
    To show us so much gentry and good will
    As to expend your time with us awhile
    For the supply and profit of our hope,
    Your visitation shall receive such thanks
    As fits a king's remembrance.
  Ros. Both your Majesties
    Might, by the sovereign power you have of us,
    Put your dread pleasures more into command
    Than to entreaty.
  Guil. But we both obey,
    And here give up ourselves, in the full bent,
    To lay our service freely at your feet,
    To be commanded.
  King. Thanks, Rosencrantz and gentle Guildenstern.
  Queen. Thanks, Guildenstern and gentle Rosencrantz.  
    And I beseech you instantly to visit
    My too much changed son.- Go, some of you,
    And bring these gentlemen where Hamlet is.
  Guil. Heavens make our presence and our practices
    Pleasant and helpful to him!
  Queen. Ay, amen!
                 Exeunt Rosencrantz and Guildenstern, [with some
                                                    Attendants].

                         Enter Polonius.

  Pol. Th' ambassadors from Norway, my good lord,
    Are joyfully return'd.
  King. Thou still hast been the father of good news.
  Pol. Have I, my lord? Assure you, my good liege,
    I hold my duty as I hold my soul,
    Both to my God and to my gracious king;
    And I do think- or else this brain of mine
    Hunts not the trail of policy so sure
    As it hath us'd to do- that I have found  
    The very cause of Hamlet's lunacy.
  King. O, speak of that! That do I long to hear.
  Pol. Give first admittance to th' ambassadors.
    My news shall be the fruit to that great feast.
  King. Thyself do grace to them, and bring them in.
                                                [Exit Polonius.]
    He tells me, my dear Gertrude, he hath found
    The head and source of all your son's distemper.
  Queen. I doubt it is no other but the main,
    His father's death and our o'erhasty marriage.
  King. Well, we shall sift him.

              Enter Polonius, Voltemand, and Cornelius.

    Welcome, my good friends.
    Say, Voltemand, what from our brother Norway?
  Volt. Most fair return of greetings and desires.
    Upon our first, he sent out to suppress
    His nephew's levies; which to him appear'd
    To be a preparation 'gainst the Polack,  
    But better look'd into, he truly found
    It was against your Highness; whereat griev'd,
    That so his sickness, age, and impotence
    Was falsely borne in hand, sends out arrests
    On Fortinbras; which he, in brief, obeys,
    Receives rebuke from Norway, and, in fine,
    Makes vow before his uncle never more
    To give th' assay of arms against your Majesty.
    Whereon old Norway, overcome with joy,
    Gives him three thousand crowns in annual fee
    And his commission to employ those soldiers,
    So levied as before, against the Polack;
    With an entreaty, herein further shown,
                                                [Gives a paper.]
    That it might please you to give quiet pass
    Through your dominions for this enterprise,
    On such regards of safety and allowance
    As therein are set down.
  King. It likes us well;
    And at our more consider'd time we'll read,  
    Answer, and think upon this business.
    Meantime we thank you for your well-took labour.
    Go to your rest; at night we'll feast together.
    Most welcome home!                       Exeunt Ambassadors.
  Pol. This business is well ended.
    My liege, and madam, to expostulate
    What majesty should be, what duty is,
    Why day is day, night is night, and time is time.
    Were nothing but to waste night, day, and time.
    Therefore, since brevity is the soul of wit,
    And tediousness the limbs and outward flourishes,
    I will be brief. Your noble son is mad.
    Mad call I it; for, to define true madness,
    What is't but to be nothing else but mad?
    But let that go.
  Queen. More matter, with less art.
  Pol. Madam, I swear I use no art at all.
    That he is mad, 'tis true: 'tis true 'tis pity;
    And pity 'tis 'tis true. A foolish figure!
    But farewell it, for I will use no art.  
    Mad let us grant him then. And now remains
    That we find out the cause of this effect-
    Or rather say, the cause of this defect,
    For this effect defective comes by cause.
    Thus it remains, and the remainder thus.
    Perpend.
    I have a daughter (have while she is mine),
    Who in her duty and obedience, mark,
    Hath given me this. Now gather, and surmise.
                                             [Reads] the letter.
    'To the celestial, and my soul's idol, the most beautified
      Ophelia,'-

    That's an ill phrase, a vile phrase; 'beautified' is a vile
      phrase.
    But you shall hear. Thus:
                                                        [Reads.]
    'In her excellent white bosom, these, &c.'
  Queen. Came this from Hamlet to her?
  Pol. Good madam, stay awhile. I will be faithful.     [Reads.]  

          'Doubt thou the stars are fire;
            Doubt that the sun doth move;
          Doubt truth to be a liar;
            But never doubt I love.
      'O dear Ophelia, I am ill at these numbers; I have not art to
    reckon my groans; but that I love thee best, O most best, believe
    it. Adieu.
      'Thine evermore, most dear lady, whilst this machine is to him,
                                                          HAMLET.'

    This, in obedience, hath my daughter shown me;
    And more above, hath his solicitings,
    As they fell out by time, by means, and place,
    All given to mine ear.
  King. But how hath she
    Receiv'd his love?
  Pol. What do you think of me?
  King. As of a man faithful and honourable.
  Pol. I would fain prove so. But what might you think,  
    When I had seen this hot love on the wing
    (As I perceiv'd it, I must tell you that,
    Before my daughter told me), what might you,
    Or my dear Majesty your queen here, think,
    If I had play'd the desk or table book,
    Or given my heart a winking, mute and dumb,
    Or look'd upon this love with idle sight?
    What might you think? No, I went round to work
    And my young mistress thus I did bespeak:
    'Lord Hamlet is a prince, out of thy star.
    This must not be.' And then I prescripts gave her,
    That she should lock herself from his resort,
    Admit no messengers, receive no tokens.
    Which done, she took the fruits of my advice,
    And he, repulsed, a short tale to make,
    Fell into a sadness, then into a fast,
    Thence to a watch, thence into a weakness,
    Thence to a lightness, and, by this declension,
    Into the madness wherein now he raves,
    And all we mourn for.  
  King. Do you think 'tis this?
  Queen. it may be, very like.
  Pol. Hath there been such a time- I would fain know that-
    That I have Positively said ''Tis so,'
    When it prov'd otherwise.?
  King. Not that I know.
  Pol. [points to his head and shoulder] Take this from this, if this
      be otherwise.
    If circumstances lead me, I will find
    Where truth is hid, though it were hid indeed
    Within the centre.
  King. How may we try it further?
  Pol. You know sometimes he walks four hours together
    Here in the lobby.
  Queen. So he does indeed.
  Pol. At such a time I'll loose my daughter to him.
    Be you and I behind an arras then.
    Mark the encounter. If he love her not,
    And he not from his reason fall'n thereon
    Let me be no assistant for a state,  
    But keep a farm and carters.
  King. We will try it.

                 Enter Hamlet, reading on a book.

  Queen. But look where sadly the poor wretch comes reading.
  Pol. Away, I do beseech you, both away
    I'll board him presently. O, give me leave.
                       Exeunt King and Queen, [with Attendants].
    How does my good Lord Hamlet?
  Ham. Well, God-a-mercy.
  Pol. Do you know me, my lord?
  Ham. Excellent well. You are a fishmonger.
  Pol. Not I, my lord.
  Ham. Then I would you were so honest a man.
  Pol. Honest, my lord?
  Ham. Ay, sir. To be honest, as this world goes, is to be one man
    pick'd out of ten thousand.
  Pol. That's very true, my lord.
  Ham. For if the sun breed maggots in a dead dog, being a god  
    kissing carrion- Have you a daughter?
  Pol. I have, my lord.
  Ham. Let her not walk i' th' sun. Conception is a blessing, but not
    as your daughter may conceive. Friend, look to't.
  Pol. [aside] How say you by that? Still harping on my daughter. Yet
    he knew me not at first. He said I was a fishmonger. He is far
    gone, far gone! And truly in my youth I suff'red much extremity
    for love- very near this. I'll speak to him again.- What do you
    read, my lord?
  Ham. Words, words, words.
  Pol. What is the matter, my lord?
  Ham. Between who?
  Pol. I mean, the matter that you read, my lord.
  Ham. Slanders, sir; for the satirical rogue says here that old men
    have grey beards; that their faces are wrinkled; their eyes
    purging thick amber and plum-tree gum; and that they have a
    plentiful lack of wit, together with most weak hams. All which,
    sir, though I most powerfully and potently believe, yet I hold it
    not honesty to have it thus set down; for you yourself, sir,
    should be old as I am if, like a crab, you could go backward.  
  Pol. [aside] Though this be madness, yet there is a method in't.-
   Will You walk out of the air, my lord?
  Ham. Into my grave?
  Pol. Indeed, that is out o' th' air. [Aside] How pregnant sometimes
    his replies are! a happiness that often madness hits on, which
    reason and sanity could not so prosperously be delivered of. I
    will leave him and suddenly contrive the means of meeting between
    him and my daughter.- My honourable lord, I will most humbly take
    my leave of you.
  Ham. You cannot, sir, take from me anything that I will more
    willingly part withal- except my life, except my life, except my
    life,

                    Enter Rosencrantz and Guildenstern.

  Pol. Fare you well, my lord.
  Ham. These tedious old fools!
  Pol. You go to seek the Lord Hamlet. There he is.
  Ros. [to Polonius] God save you, sir!
                                                Exit [Polonius].  
  Guil. My honour'd lord!
  Ros. My most dear lord!
  Ham. My excellent good friends! How dost thou, Guildenstern? Ah,
    Rosencrantz! Good lads, how do ye both?
  Ros. As the indifferent children of the earth.
  Guil. Happy in that we are not over-happy.
    On Fortune's cap we are not the very button.
  Ham. Nor the soles of her shoe?
  Ros. Neither, my lord.
  Ham. Then you live about her waist, or in the middle of her
    favours?
  Guil. Faith, her privates we.
  Ham. In the secret parts of Fortune? O! most true! she is a
    strumpet. What news ?
  Ros. None, my lord, but that the world's grown honest.
  Ham. Then is doomsday near! But your news is not true. Let me
    question more in particular. What have you, my good friends,
    deserved at the hands of Fortune that she sends you to prison
    hither?
  Guil. Prison, my lord?  
  Ham. Denmark's a prison.
  Ros. Then is the world one.
  Ham. A goodly one; in which there are many confines, wards, and
    dungeons, Denmark being one o' th' worst.
  Ros. We think not so, my lord.
  Ham. Why, then 'tis none to you; for there is nothing either good
    or bad but thinking makes it so. To me it is a prison.
  Ros. Why, then your ambition makes it one. 'Tis too narrow for your
    mind.
  Ham. O God, I could be bounded in a nutshell and count myself a
    king of infinite space, were it not that I have bad dreams.
  Guil. Which dreams indeed are ambition; for the very substance of
    the ambitious is merely the shadow of a dream.
  Ham. A dream itself is but a shadow.
  Ros. Truly, and I hold ambition of so airy and light a quality that
    it is but a shadow's shadow.
  Ham. Then are our beggars bodies, and our monarchs and outstretch'd
    heroes the beggars' shadows. Shall we to th' court? for, by my
    fay, I cannot reason.
  Both. We'll wait upon you.  
  Ham. No such matter! I will not sort you with the rest of my
    servants; for, to speak to you like an honest man, I am most
    dreadfully attended. But in the beaten way of friendship, what
    make you at Elsinore?
  Ros. To visit you, my lord; no other occasion.
  Ham. Beggar that I am, I am even poor in thanks; but I thank you;
    and sure, dear friends, my thanks are too dear a halfpenny. Were
    you not sent for? Is it your own inclining? Is it a free
    visitation? Come, deal justly with me. Come, come! Nay, speak.
  Guil. What should we say, my lord?
  Ham. Why, anything- but to th' purpose. You were sent for; and
    there is a kind of confession in your looks, which your modesties
    have not craft enough to colour. I know the good King and Queen
    have sent for you.
  Ros. To what end, my lord?
  Ham. That you must teach me. But let me conjure you by the rights
    of our fellowship, by the consonancy of our youth, by the
    obligation of our ever-preserved love, and by what more dear a
    better proposer could charge you withal, be even and direct with
    me, whether you were sent for or no.  
  Ros. [aside to Guildenstern] What say you?
  Ham. [aside] Nay then, I have an eye of you.- If you love me, hold
    not off.
  Guil. My lord, we were sent for.
  Ham. I will tell you why. So shall my anticipation prevent your
    discovery, and your secrecy to the King and Queen moult no
    feather. I have of late- but wherefore I know not- lost all my
    mirth, forgone all custom of exercises; and indeed, it goes so
    heavily with my disposition that this goodly frame, the earth,
    seems to me a sterile promontory; this most excellent canopy, the
    air, look you, this brave o'erhanging firmament, this majestical
    roof fretted with golden fire- why, it appeareth no other thing
    to me than a foul and pestilent congregation of vapours. What a
    piece of work is a man! how noble in reason! how infinite in
    faculties! in form and moving how express and admirable! in
    action how like an angel! in apprehension how like a god! the
    beauty of the world, the paragon of animals! And yet to me what
    is this quintessence of dust? Man delights not me- no, nor woman
    neither, though by your smiling you seem to say so.
  Ros. My lord, there was no such stuff in my thoughts.  
  Ham. Why did you laugh then, when I said 'Man delights not me'?
  Ros. To think, my lord, if you delight not in man, what lenten
    entertainment the players shall receive from you. We coted them
    on the way, and hither are they coming to offer you service.
  Ham. He that plays the king shall be welcome- his Majesty shall
    have tribute of me; the adventurous knight shall use his foil and
    target; the lover shall not sigh gratis; the humorous man shall
    end his part in peace; the clown shall make those laugh whose
    lungs are tickle o' th' sere; and the lady shall say her mind
    freely, or the blank verse shall halt fort. What players are
    they?
  Ros. Even those you were wont to take such delight in, the
    tragedians of the city.
  Ham. How chances it they travel? Their residence, both in
    reputation and profit, was better both ways.
  Ros. I think their inhibition comes by the means of the late
    innovation.
  Ham. Do they hold the same estimation they did when I was in the
    city? Are they so follow'd?
  Ros. No indeed are they not.  
  Ham. How comes it? Do they grow rusty?
  Ros. Nay, their endeavour keeps in the wonted pace; but there is,
    sir, an eyrie of children, little eyases, that cry out on the top
    of question and are most tyrannically clapp'd fort. These are now
    the fashion, and so berattle the common stages (so they call
    them) that many wearing rapiers are afraid of goosequills and
    dare scarce come thither.
  Ham. What, are they children? Who maintains 'em? How are they
    escoted? Will they pursue the quality no longer than they can
    sing? Will they not say afterwards, if they should grow
    themselves to common players (as it is most like, if their means
    are no better), their writers do them wrong to make them exclaim
    against their own succession.
  Ros. Faith, there has been much to do on both sides; and the nation
    holds it no sin to tarre them to controversy. There was, for a
    while, no money bid for argument unless the poet and the player
    went to cuffs in the question.
  Ham. Is't possible?
  Guil. O, there has been much throwing about of brains.
  Ham. Do the boys carry it away?  
  Ros. Ay, that they do, my lord- Hercules and his load too.
  Ham. It is not very strange; for my uncle is King of Denmark, and
    those that would make mows at him while my father lived give
    twenty, forty, fifty, a hundred ducats apiece for his picture in
    little. 'Sblood, there is something in this more than natural, if
    philosophy could find it out.

                     Flourish for the Players.

  Guil. There are the players.
  Ham. Gentlemen, you are welcome to Elsinore. Your hands, come! Th'
    appurtenance of welcome is fashion and ceremony. Let me comply
    with you in this garb, lest my extent to the players (which I
    tell you must show fairly outwards) should more appear like
    entertainment than yours. You are welcome. But my uncle-father
    and aunt-mother are deceiv'd.
  Guil. In what, my dear lord?
  Ham. I am but mad north-north-west. When the wind is southerly I
    know a hawk from a handsaw.
  
                            Enter Polonius.

  Pol. Well be with you, gentlemen!
  Ham. Hark you, Guildenstern- and you too- at each ear a hearer!
    That great baby you see there is not yet out of his swaddling
    clouts.
  Ros. Happily he's the second time come to them; for they say an old
    man is twice a child.
  Ham. I will prophesy he comes to tell me of the players. Mark it.-
   You say right, sir; a Monday morning; twas so indeed.
  Pol. My lord, I have news to tell you.
  Ham. My lord, I have news to tell you. When Roscius was an actor in
    Rome-
  Pol. The actors are come hither, my lord.
  Ham. Buzz, buzz!
  Pol. Upon my honour-
  Ham. Then came each actor on his ass-
  Pol. The best actors in the world, either for tragedy, comedy,
    history, pastoral, pastoral-comical, historical-pastoral,
    tragical-historical, tragical-comical-historical-pastoral; scene  
    individable, or poem unlimited. Seneca cannot be too heavy, nor
    Plautus too light. For the law of writ and the liberty, these are
    the only men.
  Ham. O Jephthah, judge of Israel, what a treasure hadst thou!
  Pol. What treasure had he, my lord?
  Ham. Why,

         'One fair daughter, and no more,
           The which he loved passing well.'

  Pol. [aside] Still on my daughter.
  Ham. Am I not i' th' right, old Jephthah?
  Pol. If you call me Jephthah, my lord, I have a daughter that I
    love passing well.
  Ham. Nay, that follows not.
  Pol. What follows then, my lord?
  Ham. Why,

           'As by lot, God wot,'

 and then, you know,
  
           'It came to pass, as most like it was.'

    The first row of the pious chanson will show you more; for look
    where my abridgment comes.

                     Enter four or five Players.

    You are welcome, masters; welcome, all.- I am glad to see thee
    well.- Welcome, good friends.- O, my old friend? Why, thy face is
    valanc'd since I saw thee last. Com'st' thou to' beard me in
    Denmark?- What, my young lady and mistress? By'r Lady, your
    ladyship is nearer to heaven than when I saw you last by the
    altitude of a chopine. Pray God your voice, like a piece of
    uncurrent gold, be not crack'd within the ring.- Masters, you are
    all welcome. We'll e'en to't like French falconers, fly at
    anything we see. We'll have a speech straight. Come, give us a
    taste of your quality. Come, a passionate speech.
  1. Play. What speech, my good lord?
  Ham. I heard thee speak me a speech once, but it was never acted;
    or if it was, not above once; for the play, I remember, pleas'd  
    not the million, 'twas caviary to the general; but it was (as I
    receiv'd it, and others, whose judgments in such matters cried in
    the top of mine) an excellent play, well digested in the scenes,
    set down with as much modesty as cunning. I remember one said
    there were no sallets in the lines to make the matter savoury,
    nor no matter in the phrase that might indict the author of
    affectation; but call'd it an honest method, as wholesome as
    sweet, and by very much more handsome than fine. One speech in't
    I chiefly lov'd. 'Twas AEneas' tale to Dido, and thereabout of it
    especially where he speaks of Priam's slaughter. If it live in
    your memory, begin at this line- let me see, let me see:

         'The rugged Pyrrhus, like th' Hyrcanian beast-'

    'Tis not so; it begins with Pyrrhus:

         'The rugged Pyrrhus, he whose sable arms,
         Black as his purpose, did the night resemble
         When he lay couched in the ominous horse,
         Hath now this dread and black complexion smear'd  
         With heraldry more dismal. Head to foot
         Now is be total gules, horridly trick'd
         With blood of fathers, mothers, daughters, sons,
         Bak'd and impasted with the parching streets,
         That lend a tyrannous and a damned light
         To their lord's murther. Roasted in wrath and fire,
         And thus o'ersized with coagulate gore,
         With eyes like carbuncles, the hellish Pyrrhus
         Old grandsire Priam seeks.'

    So, proceed you.
  Pol. Fore God, my lord, well spoken, with good accent and good
     discretion.

  1. Play. 'Anon he finds him,
      Striking too short at Greeks. His antique sword,
      Rebellious to his arm, lies where it falls,
      Repugnant to command. Unequal match'd,
      Pyrrhus at Priam drives, in rage strikes wide;
      But with the whiff and wind of his fell sword  
      Th' unnerved father falls. Then senseless Ilium,
      Seeming to feel this blow, with flaming top
      Stoops to his base, and with a hideous crash
      Takes prisoner Pyrrhus' ear. For lo! his sword,
      Which was declining on the milky head
      Of reverend Priam, seem'd i' th' air to stick.
      So, as a painted tyrant, Pyrrhus stood,
      And, like a neutral to his will and matter,
      Did nothing.
      But, as we often see, against some storm,
      A silence in the heavens, the rack stand still,
      The bold winds speechless, and the orb below
      As hush as death- anon the dreadful thunder
      Doth rend the region; so, after Pyrrhus' pause,
      Aroused vengeance sets him new awork;
      And never did the Cyclops' hammers fall
      On Mars's armour, forg'd for proof eterne,
      With less remorse than Pyrrhus' bleeding sword
      Now falls on Priam.
      Out, out, thou strumpet Fortune! All you gods,  
      In general synod take away her power;
      Break all the spokes and fellies from her wheel,
      And bowl the round nave down the hill of heaven,
      As low as to the fiends!

  Pol. This is too long.
  Ham. It shall to the barber's, with your beard.- Prithee say on.
    He's for a jig or a tale of bawdry, or he sleeps. Say on; come to
    Hecuba.

  1. Play. 'But who, O who, had seen the mobled queen-'

  Ham. 'The mobled queen'?
  Pol. That's good! 'Mobled queen' is good.

  1. Play. 'Run barefoot up and down, threat'ning the flames
      With bisson rheum; a clout upon that head
      Where late the diadem stood, and for a robe,
      About her lank and all o'erteemed loins,
      A blanket, in the alarm of fear caught up-  
      Who this had seen, with tongue in venom steep'd
      'Gainst Fortune's state would treason have pronounc'd.
      But if the gods themselves did see her then,
      When she saw Pyrrhus make malicious sport
      In Mincing with his sword her husband's limbs,
      The instant burst of clamour that she made
      (Unless things mortal move them not at all)
      Would have made milch the burning eyes of heaven
      And passion in the gods.'

  Pol. Look, whe'r he has not turn'd his colour, and has tears in's
    eyes. Prithee no more!
  Ham. 'Tis well. I'll have thee speak out the rest of this soon.-
    Good my lord, will you see the players well bestow'd? Do you
    hear? Let them be well us'd; for they are the abstract and brief
    chronicles of the time. After your death you were better have a
    bad epitaph than their ill report while you live.
  Pol. My lord, I will use them according to their desert.
  Ham. God's bodykins, man, much better! Use every man after his
    desert, and who should scape whipping? Use them after your own  
    honour and dignity. The less they deserve, the more merit is in
    your bounty. Take them in.
  Pol. Come, sirs.
  Ham. Follow him, friends. We'll hear a play to-morrow.
                 Exeunt Polonius and Players [except the First].
    Dost thou hear me, old friend? Can you play 'The Murther of
    Gonzago'?
  1. Play. Ay, my lord.
  Ham. We'll ha't to-morrow night. You could, for a need, study a
    speech of some dozen or sixteen lines which I would set down and
    insert in't, could you not?
  1. Play. Ay, my lord.
  Ham. Very well. Follow that lord- and look you mock him not.
                                            [Exit First Player.]
    My good friends, I'll leave you till night. You are welcome to
    Elsinore.
  Ros. Good my lord!
  Ham. Ay, so, God b' wi' ye!
                            [Exeunt Rosencrantz and Guildenstern
    Now I am alone.  
    O what a rogue and peasant slave am I!
    Is it not monstrous that this player here,
    But in a fiction, in a dream of passion,
    Could force his soul so to his own conceit
    That, from her working, all his visage wann'd,
    Tears in his eyes, distraction in's aspect,
    A broken voice, and his whole function suiting
    With forms to his conceit? And all for nothing!
    For Hecuba!
    What's Hecuba to him, or he to Hecuba,
    That he should weep for her? What would he do,
    Had he the motive and the cue for passion
    That I have? He would drown the stage with tears
    And cleave the general ear with horrid speech;
    Make mad the guilty and appal the free,
    Confound the ignorant, and amaze indeed
    The very faculties of eyes and ears.
    Yet I,
    A dull and muddy-mettled rascal, peak
    Like John-a-dreams, unpregnant of my cause,  
    And can say nothing! No, not for a king,
    Upon whose property and most dear life
    A damn'd defeat was made. Am I a coward?
    Who calls me villain? breaks my pate across?
    Plucks off my beard and blows it in my face?
    Tweaks me by th' nose? gives me the lie i' th' throat
    As deep as to the lungs? Who does me this, ha?
    'Swounds, I should take it! for it cannot be
    But I am pigeon-liver'd and lack gall
    To make oppression bitter, or ere this
    I should have fatted all the region kites
    With this slave's offal. Bloody bawdy villain!
    Remorseless, treacherous, lecherous, kindless villain!
    O, vengeance!
    Why, what an ass am I! This is most brave,
    That I, the son of a dear father murther'd,
    Prompted to my revenge by heaven and hell,
    Must (like a whore) unpack my heart with words
    And fall a-cursing like a very drab,
    A scullion!  
    Fie upon't! foh! About, my brain! Hum, I have heard
    That guilty creatures, sitting at a play,
    Have by the very cunning of the scene
    Been struck so to the soul that presently
    They have proclaim'd their malefactions;
    For murther, though it have no tongue, will speak
    With most miraculous organ, I'll have these Players
    Play something like the murther of my father
    Before mine uncle. I'll observe his looks;
    I'll tent him to the quick. If he but blench,
    I know my course. The spirit that I have seen
    May be a devil; and the devil hath power
    T' assume a pleasing shape; yea, and perhaps
    Out of my weakness and my melancholy,
    As he is very potent with such spirits,
    Abuses me to damn me. I'll have grounds
    More relative than this. The play's the thing
    Wherein I'll catch the conscience of the King.         Exit.





ACT III. Scene I.
Elsinore. A room in the Castle.

Enter King, Queen, Polonius, Ophelia, Rosencrantz, Guildenstern, and Lords.

  King. And can you by no drift of circumstance
    Get from him why he puts on this confusion,
    Grating so harshly all his days of quiet
    With turbulent and dangerous lunacy?
  Ros. He does confess he feels himself distracted,
    But from what cause he will by no means speak.
  Guil. Nor do we find him forward to be sounded,
    But with a crafty madness keeps aloof
    When we would bring him on to some confession
    Of his true state.
  Queen. Did he receive you well?
  Ros. Most like a gentleman.
  Guil. But with much forcing of his disposition.
  Ros. Niggard of question, but of our demands
    Most free in his reply.
  Queen. Did you assay him  
    To any pastime?
  Ros. Madam, it so fell out that certain players
    We o'erraught on the way. Of these we told him,
    And there did seem in him a kind of joy
    To hear of it. They are here about the court,
    And, as I think, they have already order
    This night to play before him.
  Pol. 'Tis most true;
    And he beseech'd me to entreat your Majesties
    To hear and see the matter.
  King. With all my heart, and it doth much content me
    To hear him so inclin'd.
    Good gentlemen, give him a further edge
    And drive his purpose on to these delights.
  Ros. We shall, my lord.
                            Exeunt Rosencrantz and Guildenstern.
  King. Sweet Gertrude, leave us too;
    For we have closely sent for Hamlet hither,
    That he, as 'twere by accident, may here
    Affront Ophelia.  
    Her father and myself (lawful espials)
    Will so bestow ourselves that, seeing unseen,
    We may of their encounter frankly judge
    And gather by him, as he is behav'd,
    If't be th' affliction of his love, or no,
    That thus he suffers for.
  Queen. I shall obey you;
    And for your part, Ophelia, I do wish
    That your good beauties be the happy cause
    Of Hamlet's wildness. So shall I hope your virtues
    Will bring him to his wonted way again,
    To both your honours.
  Oph. Madam, I wish it may.
                                                   [Exit Queen.]
  Pol. Ophelia, walk you here.- Gracious, so please you,
    We will bestow ourselves.- [To Ophelia] Read on this book,
    That show of such an exercise may colour
    Your loneliness.- We are oft to blame in this,
    'Tis too much prov'd, that with devotion's visage
    And pious action we do sugar o'er  
    The Devil himself.
  King. [aside] O, 'tis too true!
    How smart a lash that speech doth give my conscience!
    The harlot's cheek, beautied with plast'ring art,
    Is not more ugly to the thing that helps it
    Than is my deed to my most painted word.
    O heavy burthen!
  Pol. I hear him coming. Let's withdraw, my lord.
                                      Exeunt King and Polonius].

                           Enter Hamlet.

public:
    /**
     * Construct an exception with a message string.
     *
     * @param message  the message
     * @param filename  which source file threw the exception
     * @param rLineNumber  which line number of the source file threw the exception
     */
    Exception(std::string message, std::string filename, const unsigned rLineNumber);

    /** Get the message associated with the exception
     *
     * @return The message set when the exception was thrown.
     **/
    std::string GetMessage() const;
};

#define EXCEPTION(message) throw Exception(message, __FILE__, __LINE__)

#define NEVER_REACHED EXCEPTION("Should have been impossible to reach this line of code")

// This is to cope with NDEBUG causing variables to not be used, since they are only
// used in assert()s
#ifdef NDEBUG
#define UNUSED_OPT(var) var=var
#else
#define UNUSED_OPT(var)
#endif

// This macro is handy for calling functions like system which return non-zero on error
#define EXPECT0(cmd, arg) { \
    std::string _arg = (arg); \
    int ret = cmd(_arg.c_str()); \
    if (ret != 0) { \
        EXCEPTION("Failed to execute command: " #cmd "(" + _arg + ")"); \
    } }
// Or if you don't care about errors for some reason...
#define IGNORE_RET(cmd, arg) { \
    std::string _arg = (arg); \
    int ret = cmd(_arg.c_str()); \
    ret = ret; \
    }

#endif // _EXCEPTION_HPP_
