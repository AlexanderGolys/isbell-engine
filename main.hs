{-# OPTIONS_GHC -Wno-unrecognised-pragmas #-}
{-# HLINT ignore "Redundant bracket" #-}
{-# HLINT ignore "Use camelCase" #-}
eval :: ((a -> b), a) -> b
eval (f, x) = f x

uncurry :: (a -> b -> c) -> (a, b) -> c
uncurry f (x, y) = f x y

curry :: ((a, b) -> c) -> a -> b -> c
curry f x y = f (x, y)

splitCodom :: (a -> (b, c)) -> (a -> b, a -> c)
splitCodom f = (fst . f, snd . f)

mergeCodom :: (a -> b, a -> c) -> (a -> (b, c))
mergeCodom (f, g) x = (f x, g x)

-- Hom :: * -> * -> *

type Hom c d = c -> d
type Hom_op c d = d -> c
type Endo c = Hom c c

h_cov :: Hom a b -> (Hom x a -> Hom x b)
h_cov f g = f . g

h_x :: Hom a b -> (Hom b x -> Hom a x)
h_x f g = g . f

type R = Float
type R3 = (R, R, R)
type R4 = (R, R, R, R)
type C x = (x -> R)

type P1 = R3
type P2 = R4
type P3 = R


union_smooth :: (P3 -> C((R, R)), P1 -> C(R3), P2 -> C(R3)) -> ((P1, P2, P3) -> C(R3))



-- @since base-4.8.0.0
absurd :: Void -> a
absurd a = case a of {}

-- | If 'Void' is uninhabited then any 'Functor' that holds only
-- values of type 'Void' is holding no values.
-- It is implemented in terms of @fmap absurd@.
--
-- @since base-4.8.0.0
vacuous :: Functor f => f Void -> f a
vacuous = fmap absurd